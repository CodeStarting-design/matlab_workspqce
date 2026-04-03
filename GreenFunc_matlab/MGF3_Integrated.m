classdef MGF3_Integrated < handle
    % MGF3_Integrated - 统一的多层格林函数计算框架
    % 支持三种方法：
    %   - INTEGRATE: 直接数值积分
    %   - DCIM: 离散复镜像法
    %   - QUASISTATIC: 仅准静态项
    
    properties
        lm          % LayerManager
        smgf        % SpectralMGF
        qmgf        % QuasistaticMGF2
        dcim        % DCIM_Integrated
        f           % 频率
        s           % 设置结构
        initialized = false
        
        i_src, m_obs
        layers_set = false
    end
    
    methods
        function Initialize(obj, f, lm, s)
            % 初始化MGF引擎
            % s.method: 'INTEGRATE', 'DCIM', 'QUASISTATIC'
            % s.components: [Gxx, Gxz, Gzx, Gzz, Gphi]
            % s.extract_quasistatic: 是否提取准静态项
            % s.extract_homogeneous: 是否提取均匀项
            
            obj.f = f;
            obj.lm = lm;
            obj.s = s;
            
            % 初始化谱域格林函数
            obj.smgf = SpectralMGF();
            obj.smgf.Initialize(lm, f, s.components, ...
                isfield(s, 'extract_quasistatic') && s.extract_quasistatic, ...
                isfield(s, 'extract_singularities') && s.extract_singularities);
            
            % 初始化准静态项
            obj.qmgf = QuasistaticMGF2();
            obj.qmgf.Initialize(lm, f, s.components, ...
                isfield(s, 'extract_singularities') && s.extract_singularities);
            
            % 初始化DCIM（如果使用DCIM方法）
            if strcmpi(s.method, 'DCIM')
                obj.dcim = DCIM_Integrated();
                
                % 初始化DCIM引擎
                obj.dcim.Initialize(lm, f, s.components);
                
                % 设置DCIM参数（直接修改properties）
                obj.dcim.extract_quasistatic = ...
                    isfield(s, 'extract_quasistatic') && s.extract_quasistatic;
                obj.dcim.N1 = 100;  % strata默认
                obj.dcim.N2 = 100;
                obj.dcim.max_num_images = -1;  % 无限制
            end
            
            obj.initialized = true;
        end
        
        function SetLayers(obj, i, m, z_src, z_obs)
            % 设置源层和观测层
            % 注意：DCIM方法需要提供z_src和z_obs参数
            
            obj.i_src = i;
            obj.m_obs = m;
            obj.smgf.SetLayers(i, m);
            obj.qmgf.SetLayers(i, m);
            obj.layers_set = true;
            
            % DCIM方法需要生成复镜像
            if strcmpi(obj.s.method, 'DCIM')
                if nargin < 5
                    error('DCIM method requires z_src and z_obs in SetLayers.');
                end
                
                % 生成DCIM复镜像
                fprintf('\n生成DCIM复镜像...\n');
                obj.dcim.GenerateImages(i, m, z_src, z_obs);
                fprintf('DCIM复镜像生成完成\n\n');
            end
        end
        
        function [G_dyadic, G_phi] = ComputeMGF(obj, x, y, z, zp)
            % 计算空间域格林函数
            % 输入：直角坐标(x, y, z_obs, z_src)
            % 输出：并矢格林函数(1x9)和标量格林函数
            
            if ~obj.layers_set
                error('Layers not set. Call SetLayers first.');
            end
            
            rho = sqrt(x^2 + y^2);
            phi = atan2(y, x);
            cos_phi = cos(phi);
            sin_phi = sin(phi);
            
            % 根据方法选择计算方式
            if strcmpi(obj.s.method, 'QUASISTATIC')
                % 仅准静态项
                G_core = obj.ComputeQuasistaticOnly(rho, z, zp);
                
                % Formulation C 映射到并矢格林函数
                G_dyadic = zeros(1, 9);
                G_dyadic(1) = G_core(1);  % Gxx
                G_dyadic(3) = G_core(2) * cos_phi;  % Gxz
                G_dyadic(5) = G_core(1);  % Gyy
                G_dyadic(6) = G_core(2) * sin_phi;  % Gyz
                G_dyadic(7) = G_core(3) * cos_phi;  % Gzx
                G_dyadic(8) = G_core(3) * sin_phi;  % Gzy
                G_dyadic(9) = G_core(4);  % Gzz
                G_phi = G_core(5);
                
            elseif strcmpi(obj.s.method, 'DCIM')
                % DCIM方法（直接返回G_dyadic和G_phi）
                [G_dyadic, G_phi] = obj.ComputeMGF_DCIM(x, y, z, zp);
                
            else  % 'INTEGRATE'
                % 直接数值积分
                G_core = obj.ComputeMGF_Integration(rho, z, zp);
                
                % Formulation C 映射到并矢格林函数
                G_dyadic = zeros(1, 9);
                G_dyadic(1) = G_core(1);  % Gxx
                G_dyadic(3) = G_core(2) * cos_phi;  % Gxz
                G_dyadic(5) = G_core(1);  % Gyy
                G_dyadic(6) = G_core(2) * sin_phi;  % Gyz
                G_dyadic(7) = G_core(3) * cos_phi;  % Gzx
                G_dyadic(8) = G_core(3) * sin_phi;  % Gzy
                G_dyadic(9) = G_core(4);  % Gzz
                G_phi = G_core(5);
            end
        end
        
        function G_core = ComputeQuasistaticOnly(obj, rho, z, zp)
            % 仅计算准静态项
            K_quasi = obj.qmgf.ComputeQMGF_Spatial(z, zp, rho);
            
            % 转换到G_core格式：[Gxx, Gxz, Gzx, Gzz, Gphi]
            G_core = K_quasi;
        end
        
        function G_core = ComputeMGF_Integration(obj, rho, z, zp)
            % 直接数值积分方法（完整实现）
            
            % 准静态项
            if isfield(obj.s, 'extract_quasistatic') && obj.s.extract_quasistatic
                K_quasi = obj.qmgf.ComputeQMGF_Spatial(z, zp, rho);
            else
                K_quasi = zeros(1, 5);
            end
            
            % 设置坐标（为积分做准备）
            obj.smgf.SetSourcePoint(zp);
            obj.smgf.SetObservationPoint(z);
            
            % 积分参数
            a = 1.2 * abs(obj.lm.k_max);
            
            G_si = zeros(1, 5);
            components = [1, 2, 3, 4, 5];
            orders     = [0, 1, 1, 0, 0]; % Gxx(0), Gxz(1), Gzx(1), Gzz(0), Gphi(0)
            
            % Gxx分量
            if obj.s.components(1)
                comp = components(1);
                ord = orders(1);
                
                % 计算近场 (复路径积分)
                I_nf = SommerfeldIntegrator2.IntegrateSpectralNearField(obj.smgf, rho, comp, ord, a, false);
                
                % 计算远场 (分区外推积分)
                I_ff = SommerfeldIntegrator2.IntegrateSpectralFarField(obj.smgf, rho, comp, ord, a, false);
                
                G_si(1) = I_nf + I_ff;
            end
            
            % Gxz分量
            if obj.s.components(2)
                comp = components(2);
                ord = orders(2);
                
                % 计算近场 (复路径积分)
                I_nf = SommerfeldIntegrator2.IntegrateSpectralNearField(obj.smgf, rho, comp, ord, a, false);
                
                % 计算远场 (分区外推积分)
                I_ff = SommerfeldIntegrator2.IntegrateSpectralFarField(obj.smgf, rho, comp, ord, a, false);
                
                G_si(2) = I_nf + I_ff;
            end
            
            % Gzx分量
            if obj.s.components(3)
                comp = components(3);
                ord = orders(3);
                
                % 计算近场 (复路径积分)
                I_nf = SommerfeldIntegrator2.IntegrateSpectralNearField(obj.smgf, rho, comp, ord, a, false);
                
                % 计算远场 (分区外推积分)
                I_ff = SommerfeldIntegrator2.IntegrateSpectralFarField(obj.smgf, rho, comp, ord, a, false);
                
                G_si(3) = I_nf + I_ff;
            end
            
            % Gzz分量
            if obj.s.components(4)
                comp = components(4);
                ord = orders(4);
                
                % 计算近场 (复路径积分)
                I_nf = SommerfeldIntegrator2.IntegrateSpectralNearField(obj.smgf, rho, comp, ord, a, false);
                
                % 计算远场 (分区外推积分)
                I_ff = SommerfeldIntegrator2.IntegrateSpectralFarField(obj.smgf, rho, comp, ord, a, false);
                
                G_si(4) = I_nf + I_ff;
            end
            
            % Gphi分量
            if obj.s.components(5)
                comp = components(5);
                ord = orders(5);
                
                % 计算近场 (复路径积分)
                I_nf = SommerfeldIntegrator2.IntegrateSpectralNearField(obj.smgf, rho, comp, ord, a, false);
                
                % 计算远场 (分区外推积分)
                I_ff = SommerfeldIntegrator2.IntegrateSpectralFarField(obj.smgf, rho, comp, ord, a, false);
                
                G_si(5) = I_nf + I_ff;
            end
            
            % 结果合并        
            G_core = K_quasi + G_si;
        end
        
        function [G_dyadic, G_phi] = ComputeMGF_DCIM(obj, x, y, z, zp)
            % DCIM方法（直接调用DCIM_Integrated.ComputeSpatialMGF）

            % 调用DCIM_Integrated.ComputeSpatialMGF
            % 它已经返回完整的G_dyadic(1x9)和G_phi
            % 并且内部已经处理了准静态项
            [G_dyadic, G_phi] = obj.dcim.ComputeSpatialMGF(x, y, z, zp);
        end
    end
end
