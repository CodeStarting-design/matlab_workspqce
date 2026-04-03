classdef LayerManager < handle
    properties
        layers = []
        % 定义扁平化属性，方便调用
        eps_top, mu_top, k_top, sigma_top, isPEC_top = false
        eps_bot, mu_bot, k_bot, sigma_bot, isPEC_bot = false
        k_min, k_max
        eps  % 所有层的介电常数数组（包括半空间）
        layers_processed = false;
        
        % 原始配置暂存
        cfg_top, cfg_bot
    end
    
    methods
        function obj = LayerManager()
        end
        
        function AddLayer(obj, zmin, zmax, epsr, mur, sigma, sigmamu)
            if nargin < 7, sigmamu=0; end
            newL = struct('zmin',zmin, 'zmax',zmax, 'h',zmax-zmin, ...
                          'epsr',epsr, 'mur',mur, 'sigma',sigma, 'sigmamu',sigmamu, ...
                          'eps',0, 'mu',0, 'k',0);
            if isempty(obj.layers), obj.layers = newL; else, obj.layers(end+1) = newL; end
        end
        
        function SetHalfspaces(obj, epst, mut, sigt, epsb, mub, sigb, pect, pecb)
            obj.cfg_top = struct('eps',epst,'mu',mut,'sig',sigt,'pec',pect);
            obj.cfg_bot = struct('eps',epsb,'mu',mub,'sig',sigb,'pec',pecb);
            obj.isPEC_top = pect; obj.isPEC_bot = pecb;
        end
        
        function ProcessLayers(obj, f)
            w = 2*pi*f; eps0=8.854187817e-12; mu0=4*pi*1e-7;
            
            % Top
            obj.eps_top = eps0*obj.cfg_top.eps - 1i*obj.cfg_top.sig/w;
            obj.mu_top  = mu0*obj.cfg_top.mu;
            obj.k_top   = w * sqrt(obj.eps_top * obj.mu_top);
            
            % Bot
            obj.eps_bot = eps0*obj.cfg_bot.eps - 1i*obj.cfg_bot.sig/w;
            obj.mu_bot  = mu0*obj.cfg_bot.mu;
            obj.k_bot   = w * sqrt(obj.eps_bot * obj.mu_bot);
            
            % Layers
            for i = 1:length(obj.layers)
                L = obj.layers(i);
                L.eps = eps0*L.epsr - 1i*L.sigma/w;
                L.mu  = mu0*L.mur - 1i*L.sigmamu/w;
                L.k   = w * sqrt(L.eps * L.mu);
                obj.layers(i) = L;
            end
            
            % Stats - 只考虑层，不包含半空间（与strata C++一致）
            if ~isempty(obj.layers)
                k_vals = abs([obj.layers.k]);
            else
                k_vals = [abs(obj.k_top), abs(obj.k_bot)];
            end
            obj.k_max = max(k_vals); obj.k_min = min(k_vals);
            
            % 保存所有层的介电常数（用于DCIM）
            obj.eps = [obj.eps_top];
            if ~isempty(obj.layers)
                for i = 1:length(obj.layers)
                    obj.eps = [obj.eps, obj.layers(i).eps];
                end
            end
            obj.eps = [obj.eps, obj.eps_bot];
            
            obj.layers_processed = true;
        end
        
        % idx: -1(Top), 0..N-1(Layers), N(Bot)
        function [eps, mu, k] = getLayerParams(obj, idx)
            if idx == -1
                eps = obj.eps_top; mu = obj.mu_top; k = obj.k_top;
            elseif idx >= length(obj.layers)
                eps = obj.eps_bot; mu = obj.mu_bot; k = obj.k_bot;
            else
                L = obj.layers(idx+1); % 0-based -> 1-based
                eps = L.eps; mu = L.mu; k = L.k;
            end
        end
        
        function val = GetZmin(obj, idx)
            if idx == -1, val = obj.layers(1).zmax;
            elseif idx >= length(obj.layers), val = -inf;
            else, val = obj.layers(idx+1).zmin; end
        end
        
        function val = GetZmax(obj, idx)
            if idx == -1, val = inf;
            elseif idx >= length(obj.layers), val = obj.layers(end).zmin;
            else, val = obj.layers(idx+1).zmax; end
        end
        
        function val = GetHeight(obj, idx)
            if idx == -1 || idx >= length(obj.layers), val = inf;
            else, val = obj.layers(idx+1).h; end
        end

        function [zmin, zmax, h] = GetLayerGeom(obj, idx)
            zmin = obj.GetZmin(idx);
            zmax = obj.GetZmax(idx);
            h = obj.GetHeight(idx);
        end
        
        function idx = FindLayer(obj, z)
            if z >= obj.layers(1).zmax, idx = -1; return; end
            if z <= obj.layers(end).zmin, idx = length(obj.layers); return; end
            for i = 1:length(obj.layers)
                if z <= obj.layers(i).zmax + 1e-9 && z >= obj.layers(i).zmin - 1e-9
                    idx = i - 1; return;  % 层标号从0开始：L1=0, L2=1, L3=2, L4=3
                end
            end
            error('Layer not found for z=%g', z);
        end
    end
end