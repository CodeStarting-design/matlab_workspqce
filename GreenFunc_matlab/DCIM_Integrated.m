classdef DCIM_Integrated < handle
    % DCIM_Integrated - 基于现有GreenFunc_matlab代码的DCIM实现
    % 整合SpectralMGF和QuasistaticMGF2，参考strata的DCIM流程
    
    properties
        % 分层媒质参数
        lm          % LayerManager对象
        f           % 频率
        omega       % 角频率
        
        % 源层和观测层
        i           % 源层索引
        m           % 观测层索引
        
        % 现有引擎
        smgf        % SpectralMGF对象（已有完整Michalski-Zheng实现）
        qmgf        % QuasistaticMGF2对象（已有准静态提取实现）
        
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
        extract_quasistatic = true
        extract_singularity = false
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
            
            % 使用现有的SpectralMGF引擎
            obj.smgf = SpectralMGF();
            obj.smgf.Initialize(lm, f, obj.components, obj.extract_quasistatic, obj.extract_singularity);
            
            % 获取QuasistaticMGF2引用
            obj.qmgf = obj.smgf.qmgf;
            
            fprintf('DCIM引擎初始化完成\n');
        end
        
        function obj = GenerateImages(obj, i, m, zp, z)
            % 生成DCIM复镜像
            % 输入：
            %   i - 源层索引
            %   m - 观测层索引
            %   zp - 源点z坐标
            %   z - 观测点z坐标
            
            % 保存层索引
            obj.i = i;
            obj.m = m;
            
            fprintf('开始生成DCIM复镜像...\n');
            
            % 设置层和点
            obj.smgf.SetLayers(i, m);
            obj.smgf.SetSourcePoint(zp);
            obj.smgf.SetObservationPoint(z);
            
            % 为每个分量生成复镜像
            obj.images = struct('a', {}, 'alpha', {}, 'k', {});
            
            for comp_idx = 1:5
                if ~obj.components(comp_idx)
                    continue;
                end
                
                fprintf('  处理分量 %d...\n', comp_idx);
                
                % 生成DCIM路径（返回T02, k, t1, t2用于系数修正）
                [path1_krho, path1_kz, T02_1, k1, t1] = obj.GenerateDCIMPath_Level1(i, comp_idx);
                [path2_krho, path2_kz, T02_2, k2, t2] = obj.GenerateDCIMPath_Level2(i, comp_idx);
                
                % GPOF拟合（传入T02, k, t1, t2）
                images1 = obj.RunGPOF_Level1(path1_krho, path1_kz, T02_1, k1, t1, comp_idx, i, m, zp, z);
                images2 = obj.RunGPOF_Level2(path2_krho, path2_kz, T02_2, k2, t2, comp_idx, i, m, zp, z, images1);
                
                % 调试输出
                fprintf('    images1数量: %d\n', length(images1));
                if ~isempty(images1)
                    fprintf('    images1(1).a存在: %d, images1(1).alpha存在: %d\n', ...
                            isfield(images1(1), 'a'), isfield(images1(1), 'alpha'));
                    if isfield(images1(1), 'a')
                        fprintf('    images1(1).a = %.6e\n', images1(1).a);
                    end
                end
                
                % 合并复镜像
                % images1和images2是struct数组，每个元素有'a'和'alpha'字段
                % 需要提取所有a和alpha值
                all_a = [];
                all_alpha = [];
                
                if ~isempty(images1)
                    a1 = [images1.a];  % 提取所有a
                    alpha1 = [images1.alpha];  % 提取所有alpha
                    all_a = [all_a, a1];
                    all_alpha = [all_alpha, alpha1];
                end
                
                if ~isempty(images2)
                    a2 = [images2.a];
                    alpha2 = [images2.alpha];
                    all_a = [all_a, a2];
                    all_alpha = [all_alpha, alpha2];
                end
                
                % 使用k_min（strata使用sp.k，即k_min）
                k = obj.lm.k_min;
                
                % 保存结果
                obj.images(comp_idx).a = all_a;
                obj.images(comp_idx).alpha = all_alpha;
                obj.images(comp_idx).k = k;
                
                fprintf('    分量 %d: 生成 %d 个复镜像\n', comp_idx, length(all_a));
            end
            
            obj.is_initialized = true;
            fprintf('DCIM复镜像生成完成\n');
        end
        
        function [G_dyadic, G_phi] = ComputeSpatialMGF(obj, x, y, z, zp)
            % 计算空间格林函数（使用DCIM复镜像）
            rho = sqrt(x.^2 + y.^2);
            phi = atan2(y,x);
            cos_term = cos(phi); 
            sin_term = sin(phi);
            G = zeros(1, 5);
            
            for comp_idx = 1:5
                if ~obj.components(comp_idx)
                    continue;
                end
                
                % 获取该分量的复镜像
                a_vec = obj.images(comp_idx).a;
                alpha_vec = obj.images(comp_idx).alpha;
                k = obj.images(comp_idx).k;
                
                % 使用Sommerfeld恒等式
                G(comp_idx) = obj.ApplySommerfeldIdentity(rho, a_vec, alpha_vec, k, comp_idx);
            end
            
            % 加回准静态项
            if obj.extract_quasistatic
                % 设置准静态MGF的层索引
                obj.qmgf.SetLayers(obj.i, obj.m);
                
                % 计算空间域准静态项
                K_quasi = obj.qmgf.ComputeQMGF_Spatial(z, zp, rho);
                
                % 调试输出
                if rho > 0.009 && rho < 0.011
                    fprintf('      准静态项调试:\n');
                    fprintf('        K_quasi(1) = %.6e + %.6ei\n', real(K_quasi(1)), imag(K_quasi(1)));
                    fprintf('        G_before_quasi(1) = %.6e + %.6ei\n', real(G(1)), imag(G(1)));
                end
                
                % 加回准静态项
                G = G + K_quasi;
                
%                 % 调试输出
%                 if rho > 0.009 && rho < 0.011
%                     fprintf('        G_after_quasi(1) = %.6e + %.6ei\n', real(G(1)), imag(G(1)));
%                 end
            end

            % Formulation C 映射
            G_dyadic = zeros(1,9);
            G_dyadic(1) = G(1); % Gxx
            G_dyadic(3) = G(2) * cos_term; % Gxz
            G_dyadic(5) = G(1); % Gyy
            G_dyadic(6) = G(2) * sin_term; % Gyz
            G_dyadic(7) = G(3) * cos_term; % Gzx
            G_dyadic(8) = G(3) * sin_term; % Gzy
            G_dyadic(9) = G(4); % Gzz
            G_phi = G(5);
        end
        
        %% ========== 内部函数 ==========
        
        function [path_krho, path_kz, T02, k, t1] = GenerateDCIMPath_Level1(obj, i, comp_idx)
            % 生成DCIM路径1（远场部分）
            % 参考：strata DCIM.cpp GenerateSamplePoints_TwoLevel
            
            % 使用最小波数（strata使用k_min）
            k = obj.lm.k_min;
            
            % 计算T02（基于最大相对介电常数）
            epsr_list = [obj.lm.layers.epsr];  % 提取所有层的相对介电常数
            epsr_max = max(epsr_list);
            T02 = 5.0 * sqrt(epsr_max);
            T01 = 200.0;
            
            % 路径1：kz = -j*k*(T02 + t1), t1 ∈ [0, T01]
            t1 = linspace(0, T01, obj.N1);
            path_kz = -1j * k * (T02 + t1);
            path_krho = sqrt(k^2 - path_kz.^2);
            
            % 确保krho在第一象限
            path_krho = abs(real(path_krho)) + 1j * abs(imag(path_krho));
        end
        
        function [path_krho, path_kz, T02, k, t2] = GenerateDCIMPath_Level2(obj, i, comp_idx)
            % 生成DCIM路径2（近场部分）
            % 参考：strata DCIM.cpp GenerateSamplePoints_TwoLevel
            
            % 使用最小波数
            k = obj.lm.k_min;
            
            % 计算T02（与路径1相同）
            epsr_list = [obj.lm.layers.epsr];  % 提取所有层的相对介电常数
            epsr_max = max(epsr_list);
            T02 = 5.0 * sqrt(epsr_max);
            
            % 路径2：kz = k*(-j*t2 + (1 - t2/T02)), t2 ∈ [0.01, T02]
            t2 = linspace(1e-2, T02, obj.N2);
            path_kz = k * (-1j * t2 + (1 - t2/T02));
            path_krho = sqrt(k^2 - path_kz.^2);
            
            % 确保krho在第一象限
            path_krho = abs(real(path_krho)) + 1j * abs(imag(path_krho));
        end
        
        function images = RunGPOF_Level1(obj, path_krho, path_kz, T02, k, t1, comp_idx, i, m, zp, z)
            % 对路径1的谱域格林函数进行GPOF拟合
            % 参考：strata DCIM.cpp GenerateImages_TwoLevel
            
            % 计算谱域格林函数沿路径的值
            K_path = zeros(size(path_kz));
            
            for ii = 1:length(path_kz)
                krho = path_krho(ii);

                % 设置径向波数
                obj.smgf.SetRadialWaveNumber(krho);

                % 计算谱域格林函数
                % 注意：如果obj.smgf.extract_quasistatic=true，已经自动减去了准静态项
                obj.smgf.ComputeSpectralMGF();

                % 获取指定分量并乘以kz（关键！参考strata DCIM.cpp第804行）
                K_path(ii) = obj.smgf.K(comp_idx) * path_kz(ii);
            end

            % 注意：不需要再手动减去准静态项，因为obj.smgf已经在ComputeSpectralMGF中处理了
            
            % GPOF拟合：K(t) ~ sum a_i * exp(alpha_i * t)
            % 注意：strata使用t作为变量，而不是kz
            [alpha_t, a] = obj.GPOF_Fit_Stable(t1, K_path, obj.L1, obj.max_num_images);
            
            % 修正系数：从t空间转换到kz空间
            % 参考：strata DCIM.cpp 第406-411行
            % a_new = a * exp(-T02 * alpha_t)
            % alpha_new = alpha_t / (j*k)
            if ~isempty(a)
                a = a .* exp(-T02 * alpha_t);
                alpha = alpha_t ./ (1j * k);  % 逐元素除法
            else
                alpha = [];
            end
            
            % 保存为结构体（统一格式）
            if isempty(a)
                images = struct('a', {}, 'alpha', {});
            else
                images = struct('a', num2cell(a(:)), 'alpha', num2cell(alpha(:)));
            end
        end
        
        function images = RunGPOF_Level2(obj, path_krho, path_kz, T02, k, t2, comp_idx, i, m, zp, z, images1)
            % 对路径2进行GPOF拟合，减去Level 1的贡献
            % 参考：strata DCIM.cpp GenerateImages_TwoLevel
            
            % 计算谱域格林函数
            K_path = zeros(size(path_kz));
            
            for ii = 1:length(path_kz)
                krho = path_krho(ii);
                obj.smgf.SetRadialWaveNumber(krho);
                obj.smgf.ComputeSpectralMGF();
                K_path(ii) = obj.smgf.K(comp_idx) * path_kz(ii);
            end

            % 注意：不需要再手动减去准静态项，因为obj.smgf已经在ComputeSpectralMGF中处理了
            
            % 减去Level 1贡献
            % 参考：strata DCIM.cpp 第417-424行
            % 注意：这里使用的是修正后的a和alpha（kz空间）
            if ~isempty(images1)
                for ii = 1:length(path_kz)
                    kz = path_kz(ii);
                    for jj = 1:length(images1)
                        a = images1(jj).a;
                        alpha = images1(jj).alpha;
                        K_path(ii) = K_path(ii) - a * exp(-alpha * kz);
                    end
                end
            end
            
            % GPOF拟合剩余部分：K(t) ~ sum a_i * exp(alpha_i * t)
            [alpha_t, a] = obj.GPOF_Fit_Stable(t2, K_path, obj.L2, obj.max_num_images);
            
            % 修正系数：从t空间转换到kz空间
            % 参考：strata DCIM.cpp 第431-436行
            % 注意：Level 2的修正公式不同！
            % alpha_new = alpha_t * T02 / ((1 + j*T02) * k)
            % a_new = a * exp(k * alpha_new)
            if ~isempty(a)
                alpha = alpha_t .* T02 ./ ((1 + 1j*T02) * k);  % 逐元素运算
                a = a .* exp(k * alpha);
            else
                alpha = [];
            end
            
            % 保存为结构体
            if isempty(a)
                images = struct('a', {}, 'alpha', {});
            else
                images = struct('a', num2cell(a(:)), 'alpha', num2cell(alpha(:)));
            end
        end
        
        function [alpha, a] = GPOF_Fit_Stable(obj, t, f, L, max_images)
            % 使用现有的GPOF函数（与strata C++一致，不做额外过滤）
            % 返回：
            %   alpha - 指数（特征值）
            %   a - 系数

            if nargin < 5
                max_images = obj.max_num_images;  % 使用默认值
            end

            try
                % 调用现有GPOF函数
                % GPOF返回: [exponents, coefficients, residual]
                % 其中 y(t) = sum a_i * exp(alpha_i * t)
                
                % 与strata一致：当max_images > 0时，用max_images作为L参数
                L_use = L;
                if max_images > 0
                    L_use = max_images;
                end
                
                [alpha, a, ~] = GPOF(f(:), t(:), L_use, obj.tol_svd, obj.tol_eig);

                % 注意：strata C++ 中 RunGPOF 没有对 alpha 做任何过滤
                % 只有 SVD 截断和特征值截断（在 GPOF 函数内部完成）
                % 所以这里不做额外过滤

            catch ME
                warning(ME.identifier, 'GPOF拟合失败: %s', ME.message);
                alpha = [];
                a = [];
            end
        end
        
        function G = ApplySommerfeldIdentity(obj, rho, a_vec, alpha_vec, k, comp_idx)
            % 应用Sommerfeld恒等式
            % G(rho) = (j/2pi) * sum a_i * exp(-j*k*Rc_i) / Rc_i
            
            G = 0;
            rho_sq = rho^2;
            
            for ii = 1:length(a_vec)
                a = a_vec(ii);
                alpha = alpha_vec(ii);
                
                % 计算Rc（复距离）
                Rc = sqrt(rho_sq - alpha^2);
                
                % Sommerfeld恒等式
                if comp_idx == 1 || comp_idx == 4 || comp_idx == 5
                    % Gxx, Gzz, Gphi: 标准形式
                    G = G + a * exp(-1j * k * Rc) / Rc;
                else
                    % Gxz, Gzx: 包含导数项
                    % 参考：strata MGF.cpp 732行
                    G = G + a * rho * (1 + 1j*k*Rc) * exp(-1j*k*Rc) / (Rc^3);
                end
            end
            
            % 归一化因子
            G = G * 1j / (2 * pi);
            
            % 调试输出
            if rho > 0.009 && rho < 0.011 && comp_idx == 1
                fprintf('      ApplySommerfeldIdentity调试:\n');
                fprintf('        rho = %.4f\n', rho);
                fprintf('        k = %.6f\n', k);
                fprintf('        num_images = %d\n', length(a_vec));
                fprintf('        a(1) = %.6e, alpha(1) = %.6e + %.6ei\n', ...
                        a_vec(1), real(alpha_vec(1)), imag(alpha_vec(1)));
                fprintf('        G (before normalization) = %.6e + %.6ei\n', ...
                        real(G * 2*pi/1j), imag(G * 2*pi/1j));
                fprintf('        G (after normalization) = %.6e + %.6ei\n', ...
                        real(G), imag(G));
            end
        end
    end
end
