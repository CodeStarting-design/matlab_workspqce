classdef Layer < handle
    properties
        layerID     % 整数 ID
        zmax        % 上界面
        zmin        % 下界面
        h           % 厚度
        
        epsr        % 相对介电常数
        mur         % 相对磁导率
        sigma       % 电导率
        sigmamu     % 磁损耗 (代码里似乎有定义)
        
        % 复数参数 (Derived parameters)
        eps         % 复介电常数
        mu          % 复磁导率
        k           % 复波数
    end
    
    methods
        function obj = Layer(zmin, zmax, epsr, mur, sigma, sigmamu)
            if nargin > 0
                obj.zmin = zmin;
                obj.zmax = zmax;
                obj.epsr = epsr;
                obj.mur = mur;
                obj.sigma = sigma;
                if nargin < 6
                    obj.sigmamu = 0.0;
                else
                    obj.sigmamu = sigmamu;
                end
                
                obj.h = zmax - zmin;
                
                if obj.h <= 0
                    warning('Layer created with non-positive height: %f', obj.h);
                end
            end
        end
        
        % 更新频率相关的复数参数 (对应 ProcessLayers)
        function updateFrequencyParams(obj, omega)
            % complex epsilon = eps0 * epsr - j * sigma / omega
            obj.eps = Constants.eps0 * obj.epsr - 1i * obj.sigma / omega;
            % complex mu = mu0 * mur - j * sigmamu / omega
            obj.mu = Constants.mu0 * obj.mur - 1i * obj.sigmamu / omega;
            
            % k = omega * sqrt(mu * eps)
            obj.k = omega * sqrt(obj.mu * obj.eps);
        end
    end
end