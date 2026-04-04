classdef QuasistaticMGF2 < handle
    % QUASISTATICMGF: 准静态格林函数计算
    % 严格对应 C++ quasistatic_MGF.cpp 实现
    
    properties
        lm, f, omega
        i, m % 源层和观测层索引 
        ki, km, epsi, epsm, mui, mum
        components
        Rrs, Rls % 菲涅尔反射系数数组
        Rm, P, nu
        MVe, MVh, MIe, MIh
        
        K_spatial, K_spectral_val
        
        % 奇异性提取开关
        extract_singularities = false;
        
        % Curl 模式 (虽然您可能不用，但为了保持一致性保留)
        curl_GEM = true;
        
        % 缓存
        F % Singularity factors
        curlF
        
        % 临时变量
        curlK_spatial
        curlK_spectral
    end
    
    methods
        %% ================= Interface =================
        function Initialize(obj, lm, f, components, extract_sing)
            % components 参数在 C++ 中用于优化计算，这里暂存但不强制检查
            obj.lm = lm;
            obj.f = f;
            obj.omega = 2.0 * pi * f;
            obj.components = components;
            if nargin > 4
                obj.extract_singularities = extract_sing;
            end
            
            % 预计算所有界面的菲涅尔反射系数
            obj.ComputeFresnelCoefficients();
            
            obj.K_spatial = zeros(1, 5);
            obj.K_spectral_val = zeros(1, 5);
            obj.F = zeros(1, 5);
            
        end
        
        function SetLayers(obj, i, m)
            obj.i = i;
            obj.m = m;
            
            [obj.epsi, obj.mui, obj.ki] = obj.GetLayerParams(i);
            [obj.epsm, obj.mum, obj.km] = obj.GetLayerParams(m);
            
            if m < i
                % Observation layer is above source layer
                obj.Rm = obj.GetFresnelCoefficient_Upward(m - 1);
                obj.ComputeM_Upward();
                obj.P = 1.0;
            elseif m > i
                % Observation layer is below source layer
                obj.Rm = obj.GetFresnelCoefficient_Downward(m);
                obj.ComputeM_Downward();
                obj.P = -1.0;
            else
                obj.P = 0.0;
            end
        end
        
        function SetCurlToGEM(obj), obj.curl_GEM = true; end
        function SetCurlToGHJ(obj), obj.curl_GEM = false; end
        
        %% ================= Spectral Domain Drivers =================
        function K_spec = ComputeQMGF_Spectral(obj, z, zp, krho)
            obj.K_spectral_val = zeros(1, 5);
            
            if obj.i == obj.m
                obj.ComputeKii_Spectral(z, zp, krho);
            else
                obj.ComputeKmi_Spectral(z, zp, krho);
            end
            
            % 归一化 (对应 C++ 最后的乘法)
            mu0 = Constants.mu0; 
            eps0 = Constants.eps0;
            
            obj.K_spectral_val(1) = obj.K_spectral_val(1) * (obj.mui / mu0);
            obj.K_spectral_val(2) = obj.K_spectral_val(2) * (obj.mui / mu0);
            obj.K_spectral_val(3) = obj.K_spectral_val(3) * (obj.mui / mu0);
            obj.K_spectral_val(4) = obj.K_spectral_val(4) * (obj.mui / mu0);
            obj.K_spectral_val(5) = obj.K_spectral_val(5) * (eps0 / obj.epsi);
            
            K_spec = obj.K_spectral_val;
        end
        
        function ComputeKii_Spectral(obj, z, zp, krho)
            J = 1i;
            kzi = obj.ComputeAxialWaveNumber(obj.ki, krho);
            
            % C++ 逻辑: 按需计算，这里我们全部计算
            exp_terms_S0 = obj.ComputeExpTermsKii_kz(z, zp, kzi);
            exp_terms_S1 = obj.ComputeExpTermsKii_krho(z, zp, krho);
            
            Rr = obj.GetFresnelCoefficient_Upward(obj.i - 1);
            Rl = obj.GetFresnelCoefficient_Downward(obj.i);
            
            div_S0 = 2.0 * J * kzi;
            div_S1 = -2.0 * krho^2;
            
            % K0
            obj.K_spectral_val(1) = (exp_terms_S0(1) + Rl.Fh*exp_terms_S0(2) + Rr.Fh*exp_terms_S0(3) + Rl.Fh*Rr.Fh*(exp_terms_S0(4) + exp_terms_S0(5))) / div_S0;
            % K1
            obj.K_spectral_val(2) = ((Rl.Fe - Rl.Fh)*exp_terms_S1(2) + (-Rr.Fe + Rr.Fh)*exp_terms_S1(3) + (Rl.Fe*Rr.Fe - Rl.Fh*Rr.Fh)*(exp_terms_S1(4) - exp_terms_S1(5))) / div_S1;
            % K2
            obj.K_spectral_val(3) = ((-Rl.Fe + Rl.Fh)*exp_terms_S1(2) + (Rr.Fe - Rr.Fh)*exp_terms_S1(3) + (Rl.Fe*Rr.Fe - Rl.Fh*Rr.Fh)*(exp_terms_S1(4) - exp_terms_S1(5))) / div_S1;
            % K3
            obj.K_spectral_val(4) = (exp_terms_S0(1) + (-2.0*Rl.Fe + Rl.Fh)*exp_terms_S0(2) + (-2.0*Rr.Fe + Rr.Fh)*exp_terms_S0(3) + (2.0*Rl.Fe*Rr.Fe - Rl.Fh*Rr.Fh)*(exp_terms_S0(4) + exp_terms_S0(5))) / div_S0;
            % K4
            obj.K_spectral_val(5) = (exp_terms_S0(1) + Rl.Fe*exp_terms_S0(2) + Rr.Fe*exp_terms_S0(3) + Rl.Fe*Rr.Fe*(exp_terms_S0(4) + exp_terms_S0(5))) / div_S0;
        end
        
        function ComputeKmi_Spectral(obj, z, zp, krho)
            J = 1i;
            kzi = obj.ComputeAxialWaveNumber(obj.ki, krho);
            
            % C++ 逻辑: K1, K2 使用 krho-form; K0, K3, K4 使用 kz-form
            exp_terms_S1 = obj.ComputeExpTermsKmi_krho(z, zp, krho);
            exp_terms_S0 = obj.ComputeExpTermsKmi_kz(z, zp, kzi);
            
            divisor_S0 = 2.0 * J * kzi;
            
            Rr = obj.GetFresnelCoefficient_Upward(obj.i - 1);
            Rl = obj.GetFresnelCoefficient_Downward(obj.i);
            
            % K0
            obj.K_spectral_val(1) = (exp_terms_S0(1) + Rl.Fh*exp_terms_S0(2) + Rr.Fh*exp_terms_S0(3) + Rl.Fh*Rr.Fh*(exp_terms_S0(4) + exp_terms_S0(5)) + ...
                obj.Rm.Fh*(exp_terms_S0(6) + Rl.Fh*exp_terms_S0(7) + Rr.Fh*exp_terms_S0(8) + Rl.Fh*Rr.Fh*(exp_terms_S0(9) + exp_terms_S0(10)))) * obj.MVh / divisor_S0;
            
            % K1 (Cross)
            c1 = (obj.epsi / obj.epsm) * obj.MVe;
            c2 = (obj.mum / obj.mui) * obj.MVh;
            
            term_K1 = ((c1 - c2)*exp_terms_S1(1) ...
                + (-Rl.Fe*c1 + Rl.Fh*c2)*exp_terms_S1(2) ...
                + (-Rr.Fe*c1 + Rr.Fh*c2)*exp_terms_S1(3) ...
                + (Rl.Fe*Rr.Fe*c1 - Rl.Fh*Rr.Fh*c2)*(exp_terms_S1(4) + exp_terms_S1(5)) ...
                + (obj.Rm.Fe*c1 - obj.Rm.Fh*c2)*exp_terms_S1(6) ...
                + (-obj.Rm.Fe*Rl.Fe*c1 + obj.Rm.Fh*Rl.Fh*c2)*exp_terms_S1(7) ...
                + (-obj.Rm.Fe*Rr.Fe*c1 + obj.Rm.Fh*Rr.Fh*c2)*exp_terms_S1(8) ...
                + (obj.Rm.Fe*Rl.Fe*Rr.Fe*c1 - obj.Rm.Fh*Rl.Fh*Rr.Fh*c2)*(exp_terms_S1(9) + exp_terms_S1(10)));
            
            obj.K_spectral_val(2) = term_K1 * obj.P / (2.0 * krho^2);
            
            % K2 (Cross)
            c1_K2 = (obj.epsm * obj.mum / (obj.epsi * obj.mui)) * obj.MIe;
            c2_K2 = obj.MIh;
            
            term_K2 = ((c1_K2 - c2_K2)*exp_terms_S1(1) ...
                + (Rl.Fe*c1_K2 - Rl.Fh*c2_K2)*exp_terms_S1(2) ...
                + (Rr.Fe*c1_K2 - Rr.Fh*c2_K2)*exp_terms_S1(3) ...
                + (Rl.Fe*Rr.Fe*c1_K2 - Rl.Fh*Rr.Fh*c2_K2)*(exp_terms_S1(4) + exp_terms_S1(5)) ...
                + (-obj.Rm.Fe*c1_K2 + obj.Rm.Fh*c2_K2)*exp_terms_S1(6) ...
                + (-obj.Rm.Fe*Rl.Fe*c1_K2 + obj.Rm.Fh*Rl.Fh*c2_K2)*exp_terms_S1(7) ...
                + (-obj.Rm.Fe*Rr.Fe*c1_K2 + obj.Rm.Fh*Rr.Fh*c2_K2)*exp_terms_S1(8) ...
                + (-obj.Rm.Fe*Rl.Fe*Rr.Fe*c1_K2 + obj.Rm.Fh*Rl.Fh*Rr.Fh*c2_K2)*(exp_terms_S1(9) + exp_terms_S1(10)));
            
            obj.K_spectral_val(3) = term_K2 * obj.P / (2.0 * krho^2);
            
            % K3
            c = (1.0 + (obj.epsi * obj.mui) / (obj.epsm * obj.mum)) * obj.MIe;
            
            term_K3 = ((c - obj.MIh)*exp_terms_S0(1) ...
                 + (-Rl.Fe*c + Rl.Fh*obj.MIh)*exp_terms_S0(2) ...
                 + (-Rr.Fe*c + Rr.Fh*obj.MIh)*exp_terms_S0(3) ...
                 + (Rl.Fe*Rr.Fe*c - Rl.Fh*Rr.Fh*obj.MIh)*(exp_terms_S0(4) + exp_terms_S0(5)) ...
                 + (-obj.Rm.Fe*c + obj.Rm.Fh*obj.MIh)*exp_terms_S0(6) ...
                 + (obj.Rm.Fe*Rl.Fe*c - obj.Rm.Fh*Rl.Fh*obj.MIh)*exp_terms_S0(7) ...
                 + (obj.Rm.Fe*Rr.Fe*c - obj.Rm.Fh*Rr.Fh*obj.MIh)*exp_terms_S0(8) ...
                 + (-obj.Rm.Fe*Rl.Fe*Rr.Fe*c + obj.Rm.Fh*Rl.Fh*Rr.Fh*obj.MIh)*(exp_terms_S0(9) + exp_terms_S0(10)));
             
            obj.K_spectral_val(4) = (term_K3 / divisor_S0) * (obj.mum / obj.mui);
            
            % K4
            obj.K_spectral_val(5) = (exp_terms_S0(1) + Rl.Fe*exp_terms_S0(2) + Rr.Fe*exp_terms_S0(3) + Rl.Fe*Rr.Fe*(exp_terms_S0(4) + exp_terms_S0(5)) + ...
                obj.Rm.Fe*(exp_terms_S0(6) + Rl.Fe*exp_terms_S0(7) + Rr.Fe*exp_terms_S0(8) + Rl.Fe*Rr.Fe*(exp_terms_S0(9) + exp_terms_S0(10)))) * obj.MVe / divisor_S0;
        end
        
        %% ================= Spatial Domain Drivers =================
        function K = ComputeQMGF_Spatial(obj, z, zp, rho)
            obj.K_spatial = zeros(1, 5);
            
            if obj.i == obj.m
                obj.ComputeKii_Spatial(z, zp, rho);
            else
                obj.ComputeKmi_Spatial(z, zp, rho);
            end
            
            mu0 = Constants.mu0; 
            eps0 = Constants.eps0;
            factor = 1.0 / (4.0 * pi);
            
            % Normalize (C++ logic)
            obj.K_spatial(1) = obj.K_spatial(1) * (obj.mui / mu0) * factor;
            obj.K_spatial(2) = obj.K_spatial(2) * (obj.mui / mu0) * factor;
            obj.K_spatial(3) = obj.K_spatial(3) * (obj.mui / mu0) * factor;
            obj.K_spatial(4) = obj.K_spatial(4) * (obj.mui / mu0) * factor;
            obj.K_spatial(5) = obj.K_spatial(5) * (eps0 / obj.epsi) * factor;
            
            if obj.extract_singularities
                obj.ComputeSingularityFactors(z, zp, rho);
            end
            
            K = obj.K_spatial;
        end
        
        function ComputeKii_Spatial(obj, z, zp, rho)
            Rr = obj.GetFresnelCoefficient_Upward(obj.i - 1);
            Rl = obj.GetFresnelCoefficient_Downward(obj.i);
            if obj.extract_singularities
                obj.ComputeFii(z, zp, rho);
            end
            exp_terms_S0 = obj.ComputeExpTermsKii_SpatialS0(z, zp, rho);
            exp_terms_S1 = obj.ComputeExpTermsKii_SpatialS1(z, zp, rho);
            % K0
            obj.K_spatial(1) = exp_terms_S0(1) + Rl.Fh*exp_terms_S0(2) + Rr.Fh*exp_terms_S0(3) + Rl.Fh*Rr.Fh*(exp_terms_S0(4) + exp_terms_S0(5));
            % K1
            obj.K_spatial(2) = -1.0 * ((Rl.Fe - Rl.Fh)*exp_terms_S1(2) + (-Rr.Fe + Rr.Fh)*exp_terms_S1(3) + (Rl.Fe*Rr.Fe - Rl.Fh*Rr.Fh)*(exp_terms_S1(4) - exp_terms_S1(5)));
            % K2
            obj.K_spatial(3) = -1.0 * ((-Rl.Fe + Rl.Fh)*exp_terms_S1(2) + (Rr.Fe - Rr.Fh)*exp_terms_S1(3) + (Rl.Fe*Rr.Fe - Rl.Fh*Rr.Fh)*(exp_terms_S1(4) - exp_terms_S1(5)));
            % K3
            obj.K_spatial(4) = exp_terms_S0(1) + (-2.0*Rl.Fe + Rl.Fh)*exp_terms_S0(2) + (-2.0*Rr.Fe + Rr.Fh)*exp_terms_S0(3) + (2.0*Rl.Fe*Rr.Fe - Rl.Fh*Rr.Fh)*(exp_terms_S0(4) + exp_terms_S0(5));
            % K4
            obj.K_spatial(5) = exp_terms_S0(1) + Rl.Fe*exp_terms_S0(2) + Rr.Fe*exp_terms_S0(3) + Rl.Fe*Rr.Fe*(exp_terms_S0(4) + exp_terms_S0(5));
        end

        function ComputeFii(obj, z, zp, rho, Rr, Rl)
            % 获取菲涅尔反射系数
%             Rr = obj.GetFresnelCoefficient_Upward(obj.i - 1);
%             Rl = obj.GetFresnelCoefficient_Downward(obj.i);
        
            % 初始化默认值 (Direct Terms)
            obj.F(1) = 1.0; 
            obj.F(4) = 1.0; 
            obj.F(5) = 1.0; 
        
            obj.F(2) = 0.0; 
            obj.F(3) = 0.0; 
        
            % 检查是否在界面处重合 (Axially Coincident)
            % Case 1: 在下界面 (z == zp == zmin)
            if obj.AxiallyCoincident(z, zp) && obj.AxiallyCoincident(z, obj.lm.GetZmin(obj.i))
                obj.F(1) = obj.F(1) + Rl.Fh;
                obj.F(4) = obj.F(4) + (-2.0 * Rl.Fe + Rl.Fh);
                obj.F(5) = obj.F(5) + Rl.Fe;
        
                obj.F(2) = -1.0 * (Rl.Fe - Rl.Fh);
                obj.F(3) = -1.0 * (-Rl.Fe + Rl.Fh);
        
            % Case 2: 在上界面 (z == zp == zmax)
            elseif obj.AxiallyCoincident(z, zp) && obj.AxiallyCoincident(z, obj.lm.GetZmax(obj.i))
                obj.F(1) = obj.F(1) + Rr.Fh;
                obj.F(4) = obj.F(4) + (-2.0 * Rr.Fe + Rr.Fh);
                obj.F(5) = obj.F(5) + Rr.Fe;
        
                obj.F(2) = -1.0 * (-Rr.Fe + Rr.Fh);
                obj.F(3) = -1.0 * (Rr.Fe - Rr.Fh);
            end
        
            % 检查水平重合 (Laterally Coincident, rho -> 0)
            % 此时 J1(k*rho) -> 0，交叉项为 0
            if obj.LaterallyCoincident(rho)
                obj.F(2) = 0.0;
                obj.F(3) = 0.0;
            end
                   
            if ~obj.components(1), obj.F(1) = 0.0; end
            if ~obj.components(2), obj.F(2) = 0.0; end
            if ~obj.components(3), obj.F(3) = 0.0; end
            if ~obj.components(4), obj.F(4) = 0.0; end
            if ~obj.components(5), obj.F(5) = 0.0; end
        end
        
        function ComputeKmi_Spatial(obj, z, zp, rho)
            Rr = obj.GetFresnelCoefficient_Upward(obj.i - 1);
            Rl = obj.GetFresnelCoefficient_Downward(obj.i);
            
            if obj.extract_singularities
                obj.ComputeFmi(z, zp, rho, Rr, Rl);
            end           
            exp_terms_S0 = obj.ComputeExpTermsKmi_SpatialS0(z, zp, rho);
            exp_terms_S1 = obj.ComputeExpTermsKmi_SpatialS1(z, zp, rho);
            
            % K0
            sum0 = exp_terms_S0(1) + Rl.Fh*exp_terms_S0(2) + Rr.Fh*exp_terms_S0(3) + Rl.Fh*Rr.Fh*(exp_terms_S0(4)+exp_terms_S0(5)) + ...
                   obj.Rm.Fh*(exp_terms_S0(6) + Rl.Fh*exp_terms_S0(7) + Rr.Fh*exp_terms_S0(8) + Rl.Fh*Rr.Fh*(exp_terms_S0(9)+exp_terms_S0(10)));
            obj.K_spatial(1) = sum0 * obj.MVh;
            
            % K1
            c1 = (obj.epsi/obj.epsm)*obj.MVe; c2 = (obj.mum/obj.mui)*obj.MVh;
            sum1 = (c1 - c2)*exp_terms_S1(1) + ...
                   (-Rl.Fe*c1 + Rl.Fh*c2)*exp_terms_S1(2) + ...
                   (-Rr.Fe*c1 + Rr.Fh*c2)*exp_terms_S1(3) + ...
                   (Rl.Fe*Rr.Fe*c1 - Rl.Fh*Rr.Fh*c2)*(exp_terms_S1(4) + exp_terms_S1(5)) + ...
                   (obj.Rm.Fe*c1 - obj.Rm.Fh*c2)*exp_terms_S1(6) + ...
                   (-obj.Rm.Fe*Rl.Fe*c1 + obj.Rm.Fh*Rl.Fh*c2)*exp_terms_S1(7) + ...
                   (-obj.Rm.Fe*Rr.Fe*c1 + obj.Rm.Fh*Rr.Fh*c2)*exp_terms_S1(8) + ...
                   (obj.Rm.Fe*Rl.Fe*Rr.Fe*c1 - obj.Rm.Fh*Rl.Fh*Rr.Fh*c2)*(exp_terms_S1(9) + exp_terms_S1(10));
            obj.K_spatial(2) = sum1 * obj.P;
            
            % K2
            c1 = (obj.epsm*obj.mum/obj.epsi/obj.mui)*obj.MIe; c2 = obj.MIh;
            sum2 = (c1 - c2)*exp_terms_S1(1) + ...
                   (Rl.Fe*c1 - Rl.Fh*c2)*exp_terms_S1(2) + ...
                   (Rr.Fe*c1 - Rr.Fh*c2)*exp_terms_S1(3) + ...
                   (Rl.Fe*Rr.Fe*c1 - Rl.Fh*Rr.Fh*c2)*(exp_terms_S1(4) + exp_terms_S1(5)) + ...
                   (-obj.Rm.Fe*c1 + obj.Rm.Fh*c2)*exp_terms_S1(6) + ...
                   (-obj.Rm.Fe*Rl.Fe*c1 + obj.Rm.Fh*Rl.Fh*c2)*exp_terms_S1(7) + ...
                   (-obj.Rm.Fe*Rr.Fe*c1 + obj.Rm.Fh*Rr.Fh*c2)*exp_terms_S1(8) + ...
                   (-obj.Rm.Fe*Rl.Fe*Rr.Fe*c1 + obj.Rm.Fh*Rl.Fh*Rr.Fh*c2)*(exp_terms_S1(9) + exp_terms_S1(10));
            obj.K_spatial(3) = sum2 * obj.P;
            
            % K3
            c = (1.0 + (obj.epsi*obj.mui)/(obj.epsm*obj.mum)) * obj.MIe;
            sum3 = (c - obj.MIh)*exp_terms_S0(1) + (-Rl.Fe*c + Rl.Fh*obj.MIh)*exp_terms_S0(2) + (-Rr.Fe*c + Rr.Fh*obj.MIh)*exp_terms_S0(3) + ...
                   (Rl.Fe*Rr.Fe*c - Rl.Fh*Rr.Fh*obj.MIh)*(exp_terms_S0(4)+exp_terms_S0(5)) + ...
                   (-obj.Rm.Fe*c + obj.Rm.Fh*obj.MIh)*exp_terms_S0(6) + ...
                   (obj.Rm.Fe*Rl.Fe*c - obj.Rm.Fh*Rl.Fh*obj.MIh)*exp_terms_S0(7) + ...
                   (obj.Rm.Fe*Rr.Fe*c - obj.Rm.Fh*Rr.Fh*obj.MIh)*exp_terms_S0(8) + ...
                   (-obj.Rm.Fe*Rl.Fe*Rr.Fe*c + obj.Rm.Fh*Rl.Fh*Rr.Fh*obj.MIh)*(exp_terms_S0(9)+exp_terms_S0(10));
            obj.K_spatial(4) = sum3 * (obj.mum/obj.mui);
            
            % K4
            sum4 = exp_terms_S0(1) + Rl.Fe*exp_terms_S0(2) + Rr.Fe*exp_terms_S0(3) + Rl.Fe*Rr.Fe*(exp_terms_S0(4)+exp_terms_S0(5)) + ...
                   obj.Rm.Fe*(exp_terms_S0(6) + Rl.Fe*exp_terms_S0(7) + Rr.Fe*exp_terms_S0(8) + Rl.Fe*Rr.Fe*(exp_terms_S0(9)+exp_terms_S0(10)));
            obj.K_spatial(5) = sum4 * obj.MVe;
        end

        function ComputeFmi(obj, z, zp, rho, Rr, Rl)
            % 获取菲涅尔反射系数
%             % Upward at top of layer i-1
%             Rr = obj.GetFresnelCoefficient_Upward(obj.i - 1);
%             % Downward at bottom of layer i
%             Rl = obj.GetFresnelCoefficient_Downward(obj.i);           
            % 初始化因子 (默认值)
            obj.F(1) = obj.MVh; % (Gxx)
            
            % (Gzz)            
            term = (1.0 + (obj.epsi * obj.mui) / (obj.epsm * obj.mum)) * obj.MIe;
            obj.F(4) = term - obj.MIh; 
            
            % (Gphi)
            obj.F(5) = obj.MVe; 
            
            % Cross Terms (默认 0)
            obj.F(2) = 0.0;
            obj.F(3) = 0.0;
            
            % 检查是否轴向重合 (z == zp)
            % 仅当源和观测在界面处且重合时，奇异性行为会改变
            if obj.AxiallyCoincident(z, zp)
                % Cross Terms Update
                % F[1]
                obj.F(2) = ((obj.epsi / obj.epsm) * obj.MVe - (obj.mum / obj.mui) * obj.MVh) * obj.P;
                % F[2]
                obj.F(3) = ((obj.epsm * obj.mum / (obj.epsi * obj.mui)) * obj.MIe - obj.MIh) * obj.P;
                
                if obj.m == obj.i + 1
                    % 观测层在源层下方 (Next Layer Below)
                    obj.F(1) = (1.0 + Rl.Fh) * obj.MVh;
                    
                    % F[3]
                    term = (1.0 - Rl.Fe) * (1.0 + (obj.epsi * obj.mui) / (obj.epsm * obj.mum)) * obj.MIe;
                    obj.F(4) = term + (-1.0 + Rl.Fh) * obj.MIh;
                    
                    % F[4]
                    obj.F(5) = (1.0 + Rl.Fe) * obj.MVe;
                    
                    % Cross Terms Additions
                    obj.F(2) = obj.F(2) + (-Rl.Fe * (obj.epsi / obj.epsm) * obj.MVe + Rl.Fh * (obj.mum / obj.mui) * obj.MVh) * obj.P;
                    obj.F(3) = obj.F(3) + (Rl.Fe * (obj.epsm * obj.mum / (obj.epsi * obj.mui)) * obj.MIe - Rl.Fh * obj.MIh) * obj.P;
                    
                elseif obj.m == obj.i - 1
                    % 观测层在源层上方 (Prev Layer Above)
                    obj.F(1) = (1.0 + Rr.Fh) * obj.MVh;
                    
                    % F[3]
                    term = (1.0 - Rr.Fe) * (1.0 + (obj.epsi * obj.mui) / (obj.epsm * obj.mum)) * obj.MIe;
                    obj.F(4) = term + (-1.0 + Rr.Fh) * obj.MIh;
                    
                    % F[4]
                    obj.F(5) = (1.0 + Rr.Fe) * obj.MVe;
                    
                    % Cross Terms Additions
                    obj.F(2) = obj.F(2) + (-Rr.Fe * (obj.epsi / obj.epsm) * obj.MVe + Rr.Fh * (obj.mum / obj.mui) * obj.MVh) * obj.P;
                    obj.F(3) = obj.F(3) + (Rr.Fe * (obj.epsm * obj.mum / (obj.epsi * obj.mui)) * obj.MIe - Rr.Fh * obj.MIh) * obj.P;
                end
            end
            
            % F[3] 最终归一化
            obj.F(4) = obj.F(4) * (obj.mum / obj.mui);
            
            % 水平重合检查 (rho -> 0)
            if obj.LaterallyCoincident(rho)
                obj.F(2) = 0.0;
                obj.F(3) = 0.0;
            end
            
            if ~obj.components(1), obj.F(1) = 0.0; end
            if ~obj.components(2), obj.F(2) = 0.0; end
            if ~obj.components(3), obj.F(3) = 0.0; end
            if ~obj.components(4), obj.F(4) = 0.0; end
            if ~obj.components(5), obj.F(5) = 0.0; end            
            
        end
        
        %% ================= Helpers =================
        function exp_terms = ComputeExpTermsKii_kz(obj, z, zp, kzi)
            jkzi = 1i * kzi;
            gamma = obj.ComputeGammaTerms(z, zp);
            exp_terms = exp(-jkzi * gamma);
        end
        
        function exp_terms = ComputeExpTermsKii_krho(obj, z, zp, krho)
            gamma = obj.ComputeGammaTerms(z, zp);
            exp_terms = exp(-krho * gamma);
        end
        
        function exp_terms = ComputeExpTermsKmi_kz(obj, z, zp, kzi)
            alpha = obj.ComputeAlphaTerms(z, zp);
            exp_terms = exp(-1i * kzi * alpha);
        end
        
        function exp_terms = ComputeExpTermsKmi_krho(obj, z, zp, krho)
            alpha = obj.ComputeAlphaTerms(z, zp);
            exp_terms = exp(-krho * alpha);
        end
        
        function ksi = ComputeDistanceTermsKii(obj, z, zp, rho)
            gamma = obj.ComputeGammaTerms(z, zp);
            ksi = sqrt(rho^2 + gamma.^2);
        end
        
        function gamma = ComputeGammaTerms(obj, z, zp)
            [zmin, zmax, h] = obj.getLayerGeom(obj.i);
            gamma = [abs(z-zp), z+zp-2*zmin, 2*zmax-(z+zp), 2*h-z+zp, 2*h+z-zp];
        end
        
        function alpha = ComputeAlphaTerms(obj, z, zp)
            [zmin_i, zmax_i, h_i] = obj.getLayerGeom(obj.i);
            [zmin_m, zmax_m, ~] = obj.getLayerGeom(obj.m);
            if obj.m < obj.i, z_ = zmax_i; else, z_ = zmin_i; end
            gamma = [abs(z_-zp), z_+zp-2*zmin_i, 2*zmax_i-(z_+zp), 2*h_i-z_+zp, 2*h_i+z_-zp];
            tau = zeros(1,2);
            if obj.m < obj.i
                tau(1) = z - zmin_m; tau(2) = -z - zmin_m + 2*zmax_m;
            else
                tau(1) = zmax_m - z; tau(2) = z - 2*zmin_m + zmax_m;
            end
            alpha = zeros(1,10);
            for idx_tau = 1:2
                for idx_gamma = 1:5
                    alpha(idx_gamma + 5*(idx_tau-1)) = gamma(idx_gamma) + tau(idx_tau) + obj.nu;
                end
            end
        end
        
        % ... [Fresnel and M computations] ...
        function ComputeFresnelCoefficients(obj)
            N = length(obj.lm.layers);
            obj.Rrs = repmat(struct('Fe',0,'Fh',0), N+1, 1);
            obj.Rls = repmat(struct('Fe',0,'Fh',0), N+1, 1);
            % Upward Fresnel coefficients (Rrs)
            % C++: if (isPEC_top) Rrs[0] = {Fe=-1, Fh=-1}
            if obj.lm.isPEC_top
                obj.Rrs(1).Fh = -1;
                obj.Rrs(1).Fe = -1;
            else
                [ec, mc] = obj.GetLayerParams(0);
                [eu, mu] = obj.GetLayerParams(-1);
                obj.Rrs(1).Fh = -(mc - mu) / (mc + mu);
                obj.Rrs(1).Fe = (ec - eu) / (ec + eu);
            end
            for ii = 1:(N-1)
                [ec, mc] = obj.GetLayerParams(ii);
                [eu, mu] = obj.GetLayerParams(ii-1);
                idx = ii + 1;
                obj.Rrs(idx).Fh = -(mc - mu) / (mc + mu);
                obj.Rrs(idx).Fe = (ec - eu) / (ec + eu);
            end
            % Downward Fresnel coefficients (Rls)
            for ii = 0:(N-1)
                [ec, mc] = obj.GetLayerParams(ii);
                [el, ml] = obj.GetLayerParams(ii+1);
                idx = ii + 1;
                if ii == N-1 && obj.lm.isPEC_bot
                    obj.Rls(idx).Fh = -1;
                    obj.Rls(idx).Fe = -1;
                else
                    obj.Rls(idx).Fh = -(mc - ml) / (mc + ml);
                    obj.Rls(idx).Fe = (ec - el) / (ec + el);
                end
            end
        end
        function ComputeM_Upward(obj), obj.MVe=1; obj.MVh=1; obj.MIe=1; obj.MIh=1; obj.nu=0; for k=obj.m:(obj.i-2), Rr=obj.GetFresnelCoefficient_Upward(k); obj.MVe=obj.MVe*(1+Rr.Fe); obj.MVh=obj.MVh*(1+Rr.Fh); obj.MIe=obj.MIe*(1-Rr.Fe); obj.MIh=obj.MIh*(1-Rr.Fh); obj.nu=obj.nu+obj.lm.GetHeight(k+1); end; end
        function ComputeM_Downward(obj), obj.MVe=1; obj.MVh=1; obj.MIe=1; obj.MIh=1; obj.nu=0; for k=(obj.i+1):(obj.m-1), Rl=obj.GetFresnelCoefficient_Downward(k); obj.MVe=obj.MVe*(1+Rl.Fe); obj.MVh=obj.MVh*(1+Rl.Fh); obj.MIe=obj.MIe*(1-Rl.Fe); obj.MIh=obj.MIh*(1-Rl.Fh); obj.nu=obj.nu+obj.lm.GetHeight(k); end; end
        function f = GetFresnelCoefficient_Upward(obj, idx)
            if idx == -2
                f.Fe = 0; f.Fh = 0; return;
            end
            % 当idx=-1时，返回Rrs[1]（上半空间与第一层的Fresnel系数）
            f = obj.Rrs(idx + 2);
        end
        function f = GetFresnelCoefficient_Downward(obj, idx)
            if idx==-1 || idx==length(obj.lm.layers)
                f.Fe=0; f.Fh=0; return;
            end
            f = obj.Rls(idx + 1);
        end
        function [eps, mu, k] = GetLayerParams(obj, idx), [eps, mu, k] = obj.lm.getLayerParams(idx); end
        function [zmin, zmax, h] = getLayerGeom(obj, idx), [zmin, zmax, h] = obj.lm.GetLayerGeom(idx); end
        
        function val = GetResult_Spectral(obj, idx), val=obj.K_spectral_val(idx); end
        function val = GetResult_Spatial(obj, idx), val=obj.K_spatial(idx); end
        
        % 为了计算Kmi_SpatialS0
        function val = ComputeExpTermsKmi_SpatialS0(obj,z,zp,rho)
            alpha = obj.ComputeAlphaTerms(z, zp);
            ksi = sqrt(rho^2 + alpha.^2);
            jki = 1i * obj.ki;
            exp_terms = zeros(1,10);
            if ~obj.extract_singularities
                % 不提取奇异性：直接计算 exp(-jki*ksi)/ksi
                exp_terms = exp(-jki * ksi) ./ ksi;                
            else
                % 提取奇异性
                % 第 0 项特殊处理 (去除 1/r 奇异性)
                exp_terms(1) = obj.ComputeNonsingularHGF(ksi(1), obj.ki);
                
                if ~obj.AxiallyCoincident(z, zp)
                    % 如果 Z 轴不重合，剩下的项正常计算
                    % C++: for (int ii = 1; ii < 10; ii++)
                    exp_terms(2:10) = exp(-jki * ksi(2:10)) ./ ksi(2:10);                   
                else
                    % 如果 Z 轴重合 (z == zp)
                    if obj.m == obj.i + 1 % 观测在下一层
                        exp_terms(2) = exp_terms(1); % [1] -> (2)
                        exp_terms(3) = exp(-jki * ksi(3)) / ksi(3); % [2] -> (3)
                        
                    elseif obj.m == obj.i - 1 % 观测在上一层
                        exp_terms(2) = exp(-jki * ksi(2)) / ksi(2); % [1] -> (2)
                        exp_terms(3) = exp_terms(1); % [2] -> (3)
                    end                   
                    % 剩余项
                    % C++: for (int ii = 3; ii < 10; ii++)
                    exp_terms(4:10) = exp(-jki * ksi(4:10)) ./ ksi(4:10);
                end
            end
            val = exp_terms;
        end
        function result = ComputeNonsingularHGF(obj, r, k)
            % 计算均匀介质格林函数的非奇异部分: (exp(-jkr) - 1) / r
            % 当 r -> 0 时，通过泰勒展开避免除零错误            
            J = 1i;            
            if abs(k) < 1.0e-15
                % 波数极小 (DC 情况)，结果为 0
                result = 0.0;
                return;
            elseif abs(k * r) > 0.1
                % 距离足够远，直接使用解析公式
                result = (exp(-J * k * r) - 1.0) / r;
            else
                % 距离很近，使用泰勒级数展开
                % Series: -jk + (-jk)^2 * r / 2! + (-jk)^3 * r^2 / 3! + ...              
                threshold = 1.0e-8;          
                prev = -J * k; % 第一项 (ii=1): (-jk)^1 * r^0 / 1!
                result = prev;                
                % 从第二项开始循环 (ii=2)
                ii = 2;
                while true
                    % prev = prev * (-j*k*r) / ii
                    prev = prev * (-J * k * r) / double(ii);
                    result = result + prev;                    
                    % 收敛检查: 如果新加的项相比总结果非常小，则停止
                    if abs(prev) / abs(result) < threshold
                        break;
                    end                    
                    ii = ii + 1;
                end
            end
        end
        % 为了计算Kmi_SpatialS1
        function val = ComputeExpTermsKmi_SpatialS1(obj,z,zp,rho)
            if obj.LaterallyCoincident(rho)              
                exp_terms = zeros(1, 10); 
                val = exp_terms;
                return; % 提前退出
            end                      
            alpha = obj.ComputeAlphaTerms(z, zp);
            ksi = sqrt(rho^2 + alpha.^2);
            if ~obj.extract_singularities || ~obj.AxiallyCoincident(z, zp)
                exp_terms = (1.0 ./ rho) .* (1.0 - alpha ./ ksi);
            
            elseif obj.AxiallyCoincident(z, zp) 
                exp_terms = zeros(1, 10);              
            
                if obj.m == obj.i + 1 % C++: m == i + 1 (观测在源的下一层)                 
                    exp_terms(2) = (1.0 ./ rho) .* (1.0 - alpha(3) ./ ksi(3)); 
            
                elseif obj.m == obj.i - 1 % C++: m == i - 1 (观测在源的上一层)                  
                    exp_terms(2) = (1.0 ./ rho) .* (1.0 - alpha(2) ./ ksi(2));        
                end
                
                for ii = 4:size(exp_terms, 2) 
                    exp_terms(ii) = (1.0 ./ rho) .* (1.0 - alpha(ii) ./ ksi(ii));
                end
            end
            val = exp_terms;
        end

        % 计算Kii_SpatialS0和Kii_SpatialS1
        function exp_terms = ComputeExpTermsKii_SpatialS0(obj, z, zp, rho)
            jki = 1i * obj.ki;      
            % 1. 计算距离 ksi (array 5)
            % 假设 obj.ComputeDistanceTermsKii 已实现并返回 5 个元素
            ksi = obj.ComputeDistanceTermsKii(z, zp, rho);
            
            exp_terms = zeros(1, 5);
        
            if ~obj.extract_singularities
                % 不提取：直接计算
                exp_terms = exp(-jki * ksi) ./ ksi;
            else
                % 提取：处理奇异性
                % 第 0 项 (Direct term) 使用 Taylor 展开去除 1/r
                exp_terms(1) = obj.ComputeNonsingularHGF(ksi(1), obj.ki);
        
                % 根据 z 轴重合情况处理反射项
                % C++: if (AxiallyCoincident(z, zp) && AxiallyCoincident(z, zmin))
                if obj.AxiallyCoincident(z, zp) && obj.AxiallyCoincident(z, obj.lm.GetZmin(obj.i))
                    % 在下界面
                    exp_terms(2) = exp_terms(1); % [1] -> (2)
                    exp_terms(3) = exp(-jki * ksi(3)) / ksi(3); % [2] -> (3)
                    
                elseif obj.AxiallyCoincident(z, zp) && obj.AxiallyCoincident(z, obj.lm.GetZmax(obj.i))
                    % 在上界面
                    exp_terms(2) = exp(-jki * ksi(2)) / ksi(2); % [1] -> (2)
                    exp_terms(3) = exp_terms(1); % [2] -> (3)
                    
                else
                    % 在层中间
                    exp_terms(2) = exp(-jki * ksi(2)) / ksi(2);
                    exp_terms(3) = exp(-jki * ksi(3)) / ksi(3);
                end
        
                % 剩余项
                exp_terms(4) = exp(-jki * ksi(4)) / ksi(4);
                exp_terms(5) = exp(-jki * ksi(5)) / ksi(5);
            end
        end

        function exp_terms = ComputeExpTermsKii_SpatialS1(obj, z, zp, rho)
            % 1. 计算 ksi 和 gamma
            ksi = obj.ComputeDistanceTermsKii(z, zp, rho);
            gamma = obj.ComputeGammaTerms(z, zp);      
            exp_terms = zeros(1, 5);            
            % 如果水平距离 rho 极小，交叉项为 0 
            if obj.LaterallyCoincident(rho)
                return; % 保持全 0 返回
            end        
            if ~obj.extract_singularities
                % 不提取：正常计算                   
                for ii = 2:5
                     exp_terms(ii) = (1.0 / rho) * (1.0 - gamma(ii) / ksi(ii));
                end               
            else
                % 提取奇异性
                if obj.AxiallyCoincident(z, zp) && obj.AxiallyCoincident(z, obj.lm.GetZmin(obj.i))
                    % 下界面
                    exp_terms(2) = 0.0;
                    exp_terms(3) = (1.0 / rho) * (1.0 - gamma(3) / ksi(3));
                    
                elseif obj.AxiallyCoincident(z, zp) && obj.AxiallyCoincident(z, obj.lm.GetZmax(obj.i))
                    % 上界面
                    exp_terms(2) = (1.0 / rho) * (1.0 - gamma(2) / ksi(2));
                    exp_terms(3) = 0.0;
                    
                else
                    % 层中间
                    exp_terms(2) = (1.0 / rho) * (1.0 - gamma(2) / ksi(2));
                    exp_terms(3) = (1.0 / rho) * (1.0 - gamma(3) / ksi(3));
                end         
                % 剩余项
                exp_terms(4) = (1.0 / rho) * (1.0 - gamma(4) / ksi(4));
                exp_terms(5) = (1.0 / rho) * (1.0 - gamma(5) / ksi(5));
            end
        end

        function val = AxiallyCoincident(obj,z,zp,tol)
            if nargin < 4, tol = 1e-6; end % 默认容差    
            % 获取层叠结构的总厚度
            total_height = abs(obj.lm.layers(1).zmax - obj.lm.layers(end).zmin);           
            % 归一化容差
            tolerance = tol * total_height;            
            % 判断是否在容差范围内
            if abs(z - zp) > tolerance
                val = false;
            else
                val = true;
            end
        end
        function val = LaterallyCoincident(obj,rho,tol)
            if nargin < 3, tol = 1e-6; end % 默认容差
            % 获取层叠结构的总厚度
            total_height = abs(obj.lm.layers(1).zmax - obj.lm.layers(end).zmin);           
            % 归一化容差
            tolerance = tol * total_height;            
            % 判断 rho 是否足够小（即水平重合）
            if abs(rho) > tolerance
                val = false;
            else
                val = true;
            end           
        end
        function kz=ComputeAxialWaveNumber(obj, k, krho)
            kz=sqrt(k.^2-krho.^2); 
            % 调整kz的值，使其位于复平面的第四象限：
            % 确保波数的虚部为负，表示波在传播过程中能量衰减。
            % 确保波数的实部为正，表示波的传播方向是正向的。
            if imag(kz)>0 
                kz=complex(real(kz),-abs(imag(kz))); 
            end 
            if real(kz)<0 
                kz=complex(abs(real(kz)),imag(kz)); 
            end 
        end

        
        % 对应 C++: void QuasistaticMGF::ComputeSingularityFactors(double z, double zp, double rho)
        function ComputeSingularityFactors(obj, z, zp, rho)           
            % 初始化 F 为 0
            obj.F = zeros(1, 5);
            Rr = obj.GetFresnelCoefficient_Upward(obj.i - 1);
            Rl = obj.GetFresnelCoefficient_Downward(obj.i);
            % 根据层索引调用对应的计算函数
            if obj.i == obj.m
                obj.ComputeFii(z, zp, rho, Rr, Rl);
            else                
                obj.ComputeFmi(z, zp, rho, Rr, Rl);
            end
        
            % 归一化系数
            mu0 = Constants.mu0;
            eps0 = Constants.eps0;
            factor = 1.0 / (4.0 * pi);
        
            
            obj.F(1) = obj.F(1) * (obj.mui / mu0) * factor;
            obj.F(2) = obj.F(2) * (obj.mui / mu0) * factor;
            obj.F(3) = obj.F(3) * (obj.mui / mu0) * factor;
            obj.F(4) = obj.F(4) * (obj.mui / mu0) * factor;
            obj.F(5) = obj.F(5) * (eps0 / obj.epsi) * factor;
        
        end
    end
end