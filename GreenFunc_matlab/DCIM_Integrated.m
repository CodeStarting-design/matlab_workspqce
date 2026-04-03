classdef DCIM_Integrated < handle
    % DCIM_Integrated - 基于现有GreenFunc_matlab代码的DCIM实现
    % 整合SpectralMGF和QuasistaticMGF2，参考strata的DCIM流程
    %
    % 重要：与strata C++一致，DCIM内部的SpectralMGF使用
    % extract_quasistatic=false。DCIM直接拟合完整的谱域MGF，
    % 不需要准静态提取/加回。
    
    properties
        % 分层媒质参数
        lm          % LayerManager对象
        f           % 频率
        omega       % 角频率
        
        % 源层和观测层
        i           % 源层索引
        m           % 观测层索引
        
        % 现有引擎
        smgf        % SpectralMGF对象（extract_quasistatic=false）
        
        % DCIM参数（与strata一致）
        N1 = 100    % 路径1采样点数（strata默认100）
        N2 = 100    % 路径2采样点数（strata默认100）
        L1 = 50     % Pencil参数（路径1）
        L2 = 50     % Pencil参数（路径2）
        max_num_images = -1  % 最大复镜像数量（-1表示无限制，strata默认）
        
        % GPOF拟合容差
        tol_svd = 1e-4    % 奇异值截断容差（strata默认）
        tol_eig = 1e-16   % 特征值截断容差（strata默认）
        
        % DCIM拟合结果（每个分量一组复镜像）
        images      % 结构体数组，包含a和alpha
        
        % 设置
        components = [1, 1, 1, 1, 1];  % Gxx, Gxz, Gzx, Gzz, Gphi
        is_initialized = false
    end

    methods
        function obj = Initialize(obj, lm, f, components)
            % 初始化DCIM引擎
            obj.lm = lm;
            obj.f = f;
            obj.omega = 2*pi*f;
            
            if nargin >= 4
                obj.components = components;
            end
            
            % 关键：DCIM 内部的 SpectralMGF 使用 extract_quasistatic=false
            % 与 strata C++ testDCIM.cpp 一致：DCIM 直接拟合完整谱域 MGF
            obj.smgf = SpectralMGF();
            obj.smgf.Initialize(lm, f, obj.components, false, false);
        end
        
        function obj = GenerateImages(obj, i, m, zp, z)
            % 生成DCIM复镜像
            obj.i = i;
            obj.m = m;
            
            % 设置层和点
            obj.smgf.SetLayers(i, m);
            obj.smgf.SetSourcePoint(zp);
            obj.smgf.SetObservationPoint(z);
            
            % 为每个分量生成复镜像
            obj.images = struct('a', {}, 'alpha', {}, 'k', {});
            
            % 计算 DCIM 路径参数（共享）
            k = obj.lm.k_min;
            epsr_list = [obj.lm.layers.epsr];
            epsr_max = max(epsr_list);
            T02 = 5.0 * sqrt(epsr_max);
            T01 = 200.0;
            
            % 生成路径采样点
            t1 = linspace(0, T01, obj.N1);
            kz1 = -1j * k * (T02 + t1);
            krho1 = sqrt(k^2 - kz1.^2);
            krho1 = abs(real(krho1)) + 1j * abs(imag(krho1));
            
            t2 = linspace(1e-2, T02, obj.N2);
            kz2 = k * (-1j * t2 + (1 - t2/T02));
            krho2 = sqrt(k^2 - kz2.^2);
            krho2 = abs(real(krho2)) + 1j * abs(imag(krho2));
            
            % 采样所有分量在路径1上
            K_path1 = zeros(5, obj.N1);
            for mm = 1:obj.N1
                obj.smgf.SetRadialWaveNumber(krho1(mm));
                obj.smgf.ComputeSpectralMGF();
                for ii = 1:5
                    if obj.components(ii)
                        K_path1(ii, mm) = obj.smgf.K(ii) * kz1(mm);
                    end
                end
            end
            
            % 采样所有分量在路径2上
            K_path2 = zeros(5, obj.N2);
            for mm = 1:obj.N2
                obj.smgf.SetRadialWaveNumber(krho2(mm));
                obj.smgf.ComputeSpectralMGF();
                for ii = 1:5
                    if obj.components(ii)
                        K_path2(ii, mm) = obj.smgf.K(ii) * kz2(mm);
                    end
                end
            end
            
            % dt
            dt1 = t1(2) - t1(1);
            dt2 = t2(2) - t2(1);
            
            % 为每个分量生成复镜像
            for comp_idx = 1:5
                if ~obj.components(comp_idx)
                    obj.images(comp_idx).a = [];
                    obj.images(comp_idx).alpha = [];
                    obj.images(comp_idx).k = k;
                    continue;
                end
                
                % === Level 1: GPOF ===
                [alpha_t, a_t, ~] = GPOF(K_path1(comp_idx,:).', t1(:), obj.L1, obj.tol_svd, obj.tol_eig);
                
                % t->kz 转换 (Level 1)
                % C++: a[jj] = a[jj]*exp(-T02*alpha[jj]);
                %      alpha[jj] = alpha[jj]/(J*k);
                if ~isempty(a_t)
                    a1 = a_t .* exp(-T02 * alpha_t);
                    alpha1 = alpha_t ./ (1j * k);
                else
                    a1 = [];
                    alpha1 = [];
                end
                
                % 从路径2数据中减去 Level 1 贡献
                K2_mod = K_path2(comp_idx, :);
                for mm = 1:obj.N2
                    for jj = 1:length(a1)
                        K2_mod(mm) = K2_mod(mm) - a1(jj) * exp(-alpha1(jj) * kz2(mm));
                    end
                end
                
                % === Level 2: GPOF ===
                [alpha_t2, a_t2, ~] = GPOF(K2_mod(:), t2(:), obj.L2, obj.tol_svd, obj.tol_eig);
                
                % t->kz 转换 (Level 2)
                % C++: alpha[jj] = alpha[jj]*T02/((1.0+J*T02)*k);
                %      a[jj] = a[jj]*exp(k*alpha[jj]);
                if ~isempty(a_t2)
                    alpha2 = alpha_t2 .* T02 ./ ((1 + 1j*T02) * k);
                    a2 = a_t2 .* exp(k * alpha2);
                else
                    a2 = [];
                    alpha2 = [];
                end
                
                % 合并 Level 1 和 Level 2
                all_a = [a1(:); a2(:)];
                all_alpha = [alpha1(:); alpha2(:)];
                
                % 保存
                obj.images(comp_idx).a = all_a.';
                obj.images(comp_idx).alpha = all_alpha.';
                obj.images(comp_idx).k = k;
            end
            
            obj.is_initialized = true;
        end
        
        function [G_dyadic, G_phi] = ComputeSpatialMGF(obj, x, y, z, zp)
            % 计算空间格林函数（使用DCIM复镜像 + Sommerfeld恒等式）
            % 注意：DCIM 直接拟合完整谱域 MGF，不需要加回准静态项
            rho = sqrt(x.^2 + y.^2);
            phi = atan2(y, x);
            cos_term = cos(phi); 
            sin_term = sin(phi);
            G = zeros(1, 5);
            
            for comp_idx = 1:5
                if ~obj.components(comp_idx)
                    continue;
                end
                
                a_vec = obj.images(comp_idx).a;
                alpha_vec = obj.images(comp_idx).alpha;
                k = obj.images(comp_idx).k;
                
                if isempty(a_vec)
                    continue;
                end
                
                % Sommerfeld 恒等式
                rho_sq = rho^2;
                G_val = 0;
                for ii = 1:length(a_vec)
                    Rc = sqrt(rho_sq - alpha_vec(ii)^2);
                    if comp_idx == 1 || comp_idx == 4 || comp_idx == 5
                        % S0: Gxx, Gzz, Gphi
                        G_val = G_val + a_vec(ii) * exp(-1j * k * Rc) / Rc;
                    else
                        % S1: Gxz, Gzx
                        G_val = G_val + rho * a_vec(ii) * (1 + 1j*k*Rc) * exp(-1j*k*Rc) / (Rc^3);
                    end
                end
                G(comp_idx) = G_val * 1j / (2 * pi);
            end

            % Formulation C 映射
            G_dyadic = zeros(1, 9);
            G_dyadic(1) = G(1);              % Gxx
            G_dyadic(3) = G(2) * cos_term;   % Gxz
            G_dyadic(5) = G(1);              % Gyy
            G_dyadic(6) = G(2) * sin_term;   % Gyz
            G_dyadic(7) = G(3) * cos_term;   % Gzx
            G_dyadic(8) = G(3) * sin_term;   % Gzy
            G_dyadic(9) = G(4);              % Gzz
            G_phi = G(5);
        end
    end
end
