classdef SpectralMGF < handle
    properties
        lm, f, omega, qmgf
        components, extract_quasistatic, extract_singularities, curl_GEM=false
        initialized=false, layers_set=false, src_point_set=false, obs_point_set=false, precomputations_done=false
        i, m, z, zp, krho, krhosq
        ki, km, epsi, epsm, mui, mum
        kz0, kzi, kzm, jkzi, jkzm
        Z_top, Z_bot, Zi, Zm, Zdi, Gri, Gli, Gm, De, Dh
        MVe, MVh, MIe, MIh
        GVe, GVh, GIe, GIh 
        WVe, WVh, WIe, WIh 
        TVe, TVh, TIe, TIh
        exp_terms, trf_terms, cos_term, sin_term
        exp_terms_computed=false, trf_terms_computed=false, src_terms_computed=false, obs_terms_computed=false
        K, curlK
    end
    
    methods
        function Initialize(obj, lm, f, components, extract_quasi, extract_singu)
            obj.lm = lm; obj.f = f; obj.omega = 2*pi*f;
            if nargin < 4, components = true(1,5); end
            if nargin < 5, extract_quasi = false; end
            if nargin < 6, extract_singu = false; end
            obj.components = components;
            obj.extract_quasistatic = extract_quasi;
            obj.extract_singularities = extract_singu;
            obj.initialized = true; 
            obj.ResetState();
            if extract_quasi
                obj.qmgf = QuasistaticMGF2(); 
                obj.qmgf.Initialize(lm,f,obj.components,obj.extract_singularities);
            end
        end
        function ResetState(obj) 
            obj.layers_set=false; 
            obj.src_point_set=false; 
            obj.obs_point_set=false; 
            obj.precomputations_done=false; 
        end
        function SetLayers(obj, i, m) 
            obj.i=i; 
            obj.m=m; 
            [obj.epsi, obj.mui, obj.ki] = obj.GetLayerParams(i); 
            [obj.epsm, obj.mum, obj.km] = obj.GetLayerParams(m); 
            if (obj.extract_quasistatic)
		        obj.qmgf.SetLayers(i, m);
            end
            obj.layers_set=true; 
            obj.ResetComputedFlags(); 
        end
        function SetSourcePoint(obj, zp) 
            obj.zp=zp; 
            obj.src_point_set=true; 
            if obj.precomputations_done && obj.i~=obj.m 
                obj.src_terms_computed=false; 
            end
            obj.exp_terms_computed=false; 
        end
        function SetObservationPoint(obj, z) 
            obj.z=z; 
            obj.obs_point_set=true; 
            obj.obs_terms_computed=false; 
            obj.exp_terms_computed=false; 
            obj.trf_terms_computed=false; 
        end
        
        function SetRadialWaveNumber(obj, krho)
            if ~obj.layers_set, error('SetLayers first'); end
            obj.krho=krho; 
            obj.krhosq=krho^2;
            obj.kz0 = obj.ComputeAxialWaveNumber(obj.lm.k_min, krho);
            obj.kzi = obj.ComputeAxialWaveNumber(obj.ki, krho);
            obj.kzm = obj.ComputeAxialWaveNumber(obj.km, krho);
            obj.jkzi = 1i*obj.kzi; 
            obj.jkzm = 1i*obj.kzm;
            
            obj.ComputeHalfspaceImpedance(krho);
            obj.Zi = obj.ComputeLayerImpedance(krho, obj.i, obj.kzi);
            obj.Zm = obj.ComputeLayerImpedance(krho, obj.m, obj.kzm);
            
            % Fix: removed ComputeAllReflections, use recursive directly
            obj.Gri = obj.ComputeReflectionCoefficient_Upward(krho, obj.i-1);
            obj.Gli = obj.ComputeReflectionCoefficient_Downward(krho, obj.i);
            
            h_i = obj.lm.GetHeight(obj.i);
            if obj.i==-1 || obj.i==length(obj.lm.layers)
                exp_term=0; 
            else
                exp_term=exp(-obj.jkzi*2.0*h_i); 
            end
            obj.De = 1.0 - obj.Gli.Ge.*obj.Gri.Ge.*exp_term;
            obj.Dh = 1.0 - obj.Gli.Gh.*obj.Gri.Gh.*exp_term;
            
            if obj.m < obj.i
                obj.Zdi = obj.ComputeInterfaceImpedance_Upward(krho, obj.i-1);
                obj.Gm  = obj.ComputeReflectionCoefficient_Upward(krho, obj.m-1);
                [obj.MVe, obj.MVh, obj.MIe, obj.MIh] = obj.ComputeM_Upward(krho);
            elseif obj.m > obj.i
                obj.Zdi = obj.ComputeInterfaceImpedance_Downward(krho, obj.i);
                obj.Gm  = obj.ComputeReflectionCoefficient_Downward(krho, obj.m);
                [obj.MVe, obj.MVh, obj.MIe, obj.MIh] = obj.ComputeM_Downward(krho);
            end
            obj.precomputations_done=true; 
            obj.ResetComputedFlags();
        end
        
        function ComputeSpectralMGF(obj)
            obj.K = zeros(1,5);
            if obj.i == obj.m, obj.ComputeKii(); else, obj.ComputeKmi(); end
            J=1i; 
            mu0=Constants.mu0; 
            eps0=Constants.eps0; 
            w=obj.omega;

            obj.K(1) = obj.K(1) ./ (J*w*mu0);
            obj.K(2) = obj.K(2) .* (-obj.mui ./ (J*w*obj.epsm*mu0));
            obj.K(3) = obj.K(3) ./ (-J*w*mu0);
            obj.K(4) = obj.K(4) .* (obj.mum ./ (J * w * obj.epsi * mu0));
            obj.K(5) = obj.K(5) .* (J*w*eps0);
            
            % 这里还没调试
            if obj.extract_quasistatic
                K_quasi = obj.qmgf.ComputeQMGF_Spectral(obj.z, obj.zp, obj.krho);
                obj.K = obj.K - K_quasi;
            end

        end

        function val = GetResult(obj, idx), val = obj.K(idx); end
        
        function ComputeKii(obj)
            if obj.components(1)||obj.components(5) 
                [obj.GVe,obj.GVh]=obj.ComputeGVii(obj.z,obj.zp); 
            end
            if obj.components(1) 
                obj.K(1) = obj.GVh; 
            end
            if obj.components(2) 
                [obj.WVe,obj.WVh]=obj.ComputeWVii(obj.z,obj.zp); 
                obj.K(2)=(obj.WVe-(obj.km^2/obj.kzm^2)*obj.WVh)./obj.krhosq; 
            end
            if obj.components(3) 
                [obj.WIe,obj.WIh]=obj.ComputeWIii(obj.z,obj.zp); 
                obj.K(3)=((obj.km^2/obj.kzm^2)*obj.WIe-obj.WIh)./obj.krhosq; 
            end
            if obj.components(4) 
                [obj.GIe,obj.GIh]=obj.ComputeGIii(obj.z,obj.zp);
                obj.K(4)=obj.GIe-(obj.ki^2/obj.krhosq)*((obj.kzm^2/obj.km^2)*obj.GIe-obj.GIh);
            end
            if obj.components(5) 
                obj.K(5)=(obj.GVe-obj.GVh)./obj.krhosq; 
            end
        end
        
        function ComputeKmi(obj)                      
            if obj.m < obj.i, z_ref=obj.lm.GetZmax(obj.i); I_sign=1.0;
            elseif obj.m > obj.i, z_ref=obj.lm.GetZmin(obj.i); I_sign=-1.0; end
            if ~obj.src_terms_computed
                if obj.components(1) || obj.components(3) || obj.components(5)
                    [obj.GVe,obj.GVh]=obj.ComputeGVii(z_ref,obj.zp); 
                end
                if obj.components(2) || obj.components(4)
                    [obj.GIe,obj.GIh]=obj.ComputeGIii(z_ref,obj.zp); 
                end
                obj.src_terms_computed=true; 
            end
            if ~obj.obs_terms_computed 
                if obj.components(1) || obj.components(2) || obj.components(5)
                    obj.ComputeTV(obj.z); 
                end
                if obj.components(3) || obj.components(4)
                    obj.ComputeTI(obj.z); 
                end
                obj.obs_terms_computed=true; 
            end
            if obj.components(1) 
                obj.K(1)=obj.GVh.*obj.TVh; 
            end
            if obj.components(2) 
                IVe_mi=I_sign*obj.Zdi.Ze.*obj.GIe.*obj.TVe; 
                IVh_mi=I_sign*obj.Zdi.Zh.*obj.GIh.*obj.TVh; 
                WVe_mi=(-1i*obj.omega*obj.epsm).*IVe_mi; 
                WVh_mi=(-1i*obj.kzm.^2/(obj.omega*obj.mum)).*IVh_mi; 
                obj.K(2)=(WVe_mi-(obj.km^2./obj.kzm.^2).*WVh_mi)./obj.krhosq; 
            end
            if obj.components(3) 
                IIe_mi=I_sign*obj.Zdi.Ye.*obj.GVe.*obj.TIe; 
                IIh_mi=I_sign*obj.Zdi.Yh.*obj.GVh.*obj.TIh; 
                WIe_mi = (-1i * obj.kzm ./ obj.Zm.Ye) .* IIe_mi;
                WIh_mi = (-1i * obj.kzm ./ obj.Zm.Yh) .* IIh_mi;
                obj.K(3)=((obj.km^2./obj.kzm.^2).*WIe_mi-WIh_mi)./obj.krhosq; 
            end
            if obj.components(4) 
                GIe_mi=obj.GIe.*obj.TIe; 
                GIh_mi=obj.GIh.*obj.TIh; 
                obj.K(4)=GIe_mi-(obj.ki^2./obj.krhosq).*((obj.kzm.^2./obj.km^2).*GIe_mi-GIh_mi); 
               
            end
            if obj.components(5) 
                obj.K(5)=(obj.GVe.*obj.TVe-obj.GVh.*obj.TVh)./obj.krhosq; 
            end
            % 调试用
%             obj.LogVariables();
        end
        
        function [gve, gvh] = ComputeGVii(obj, z_val, zp_val)
            if ~obj.exp_terms_computed, obj.ComputeExpTermsGii(z_val, zp_val); end
            termE=(1./obj.De).*(obj.Gli.Ge.*obj.exp_terms(2)+obj.Gri.Ge.*obj.exp_terms(3)+2*obj.Gli.Ge.*obj.Gri.Ge.*obj.cos_term);
            gve=(obj.Zi.Ze/2).*(obj.exp_terms(1)+termE);
            termH=(1./obj.Dh).*(obj.Gli.Gh.*obj.exp_terms(2)+obj.Gri.Gh.*obj.exp_terms(3)+2*obj.Gli.Gh.*obj.Gri.Gh.*obj.cos_term);
            gvh=(obj.Zi.Zh/2).*(obj.exp_terms(1)+termH);
        end
        function [gie, gih] = ComputeGIii(obj, z_val, zp_val)
            if ~obj.exp_terms_computed, obj.ComputeExpTermsGii(z_val, zp_val); end
            termE=(1./obj.De).*(-obj.Gli.Ge.*obj.exp_terms(2)-obj.Gri.Ge.*obj.exp_terms(3)+2*obj.Gli.Ge.*obj.Gri.Ge.*obj.cos_term);
            gie=(obj.Zi.Ye/2).*(obj.exp_terms(1)+termE);
            termH=(1./obj.Dh).*(-obj.Gli.Gh.*obj.exp_terms(2)-obj.Gri.Gh.*obj.exp_terms(3)+2*obj.Gli.Gh.*obj.Gri.Gh.*obj.cos_term);
            gih=(obj.Zi.Yh/2).*(obj.exp_terms(1)+termH);
        end
        function [wve, wvh] = ComputeWVii(obj, z_val, zp_val)
            if ~obj.exp_terms_computed, obj.ComputeExpTermsGii(z_val, zp_val); end; J=1i;
            wve=(1i*obj.omega*obj.epsi/2./obj.De).*(obj.Gli.Ge.*obj.exp_terms(2)-obj.Gri.Ge.*obj.exp_terms(3)+2*J*obj.Gli.Ge.*obj.Gri.Ge.*obj.sin_term);
            wvh=(1i*obj.kzi.^2/(2*obj.omega*obj.mui)./obj.Dh).*(obj.Gli.Gh.*obj.exp_terms(2)-obj.Gri.Gh.*obj.exp_terms(3)+2*J*obj.Gli.Gh.*obj.Gri.Gh.*obj.sin_term);
        end
        function [wie, wih] = ComputeWIii(obj, z_val, zp_val)
            if ~obj.exp_terms_computed, obj.ComputeExpTermsGii(z_val, zp_val); end; J=1i;
            wie=(1i*obj.kzi.^2/(2*obj.omega*obj.epsi)./obj.De).*(-obj.Gli.Ge.*obj.exp_terms(2)+obj.Gri.Ge.*obj.exp_terms(3)+2*J*obj.Gli.Ge.*obj.Gri.Ge.*obj.sin_term);
            wih=(1i*obj.omega*obj.mui/2./obj.Dh).*(-obj.Gli.Gh.*obj.exp_terms(2)+obj.Gri.Gh.*obj.exp_terms(3)+2*J*obj.Gli.Gh.*obj.Gri.Gh.*obj.sin_term);
        end
        function ComputeExpTermsGii(obj, z_val, zp_val)
            zmin=obj.lm.GetZmin(obj.i); zmax=obj.lm.GetZmax(obj.i); h=obj.lm.GetHeight(obj.i);
            obj.exp_terms=zeros(1,3); v1=exp(-obj.jkzi*abs(z_val-zp_val)); v2=exp(-obj.jkzi*(z_val+zp_val-2.0*zmin)); v3=exp(-obj.jkzi*(2*zmax-(z_val+zp_val)));
            obj.exp_terms(1)=v1(1); obj.exp_terms(2)=v2(1); obj.exp_terms(3)=v3(1);
            if isinf(h), obj.cos_term=0; obj.sin_term=0; else, ta=exp(obj.jkzi*(-2*h+z_val-zp_val)); tb=exp(-obj.jkzi*(2*h+z_val-zp_val)); obj.cos_term=(ta(1)+tb(1))/2; obj.sin_term=(ta(1)-tb(1))/2i; end
            if obj.i==-1, obj.exp_terms(3)=0; elseif obj.i==length(obj.lm.layers), obj.exp_terms(2)=0; end
            obj.exp_terms_computed=true;
        end
        function ComputeTV(obj, z_val)
            if ~obj.trf_terms_computed, obj.ComputeTrfTerms(z_val); end
            obj.TVe=(obj.trf_terms(1)+obj.Gm.Ge.*obj.trf_terms(2)).*obj.MVe;
            obj.TVh=(obj.trf_terms(1)+obj.Gm.Gh.*obj.trf_terms(2)).*obj.MVh;
        end
        function ComputeTI(obj, z_val)
            if ~obj.trf_terms_computed
                obj.ComputeTrfTerms(z_val); 
            end
            obj.TIe=(obj.trf_terms(1)-obj.Gm.Ge.*obj.trf_terms(2)).*obj.MIe;
            obj.TIh=(obj.trf_terms(1)-obj.Gm.Gh.*obj.trf_terms(2)).*obj.MIh;
        end
        function ComputeTrfTerms(obj, z_val)
            zmin=obj.lm.GetZmin(obj.m); zmax=obj.lm.GetZmax(obj.m);
            obj.trf_terms=zeros(1,2);
            if obj.m<obj.i, v1=exp(-obj.jkzm*(z_val-zmin)); v2=exp(-obj.jkzm*(2*zmax-z_val-zmin));
            elseif obj.m>obj.i, v1=exp(-obj.jkzm*(zmax-z_val)); v2=exp(-obj.jkzm*(z_val-2*zmin+zmax)); end
            obj.trf_terms(1)=v1(1); obj.trf_terms(2)=v2(1); obj.trf_terms_computed=true;
        end
        
        function ComputeHalfspaceImpedance(obj, krho)
            kz=obj.ComputeAxialWaveNumber(obj.lm.k_top,krho); 
            % Fix: isPEC_top
            if obj.lm.isPEC_top
                obj.Z_top=struct('Ze',0,'Ye',inf,'Zh',0,'Yh',inf); 
            else 
                obj.Z_top.Ze=kz./(obj.omega*obj.lm.eps_top); 
                obj.Z_top.Ye=1./obj.Z_top.Ze; 
                obj.Z_top.Zh=(obj.omega*obj.lm.mu_top)./kz; 
                obj.Z_top.Yh=1./obj.Z_top.Zh; 
            end
            kz=obj.ComputeAxialWaveNumber(obj.lm.k_bot,krho); 
            % Fix: isPEC_bot
            if obj.lm.isPEC_bot
                obj.Z_bot=struct('Ze',0,'Ye',inf,'Zh',0,'Yh',inf); 
            else 
                obj.Z_bot.Ze=kz./(obj.omega*obj.lm.eps_bot); 
                obj.Z_bot.Ye=1./obj.Z_bot.Ze; 
                obj.Z_bot.Zh=(obj.omega*obj.lm.mu_bot)./kz; 
                obj.Z_bot.Yh=1./obj.Z_bot.Zh; 
            end
        end
        function Z=ComputeLayerImpedance(obj, krho, idx, kz)
            if idx==-1
                obj.ComputeHalfspaceImpedance(krho);
                Z=obj.Z_top; 
                return; 
            end
            if idx==length(obj.lm.layers) 
                obj.ComputeHalfspaceImpedance(krho);
                Z=obj.Z_bot; 
                return; 
            end
            [eps,mu,k]=obj.GetLayerParams(idx); 
            if nargin<4||isempty(kz)
                kz=obj.ComputeAxialWaveNumber(k,krho);             
            end
            Z.Ze=kz./(obj.omega*eps); 
            Z.Ye=1./Z.Ze; 
            Z.Zh=(obj.omega*mu)./kz; 
            Z.Yh=1./Z.Zh;
        end
        function Zr=ComputeInterfaceImpedance_Upward(obj, krho, idx)
            if idx==-1, Zr=obj.Z_top; return; end
            % Fix: isPEC_top
            if idx==0 && obj.lm.isPEC_top, [~,~,k]=obj.GetLayerParams(0); kz=obj.ComputeAxialWaveNumber(k,krho); t=tan(kz*obj.lm.GetHeight(0)); Z=obj.ComputeLayerImpedance(krho,0,kz); J=1i; Zr.Ze=Z.Ze.*J.*t; Zr.Zh=Z.Zh.*J.*t; Zr.Ye=Z.Ye./(J.*t); Zr.Yh=Z.Yh./(J.*t); return; end
            [~,~,k]=obj.GetLayerParams(idx); kz=obj.ComputeAxialWaveNumber(k,krho); t=tan(kz*obj.lm.GetHeight(idx)); Z=obj.ComputeLayerImpedance(krho,idx,kz); Zr_next=obj.ComputeInterfaceImpedance_Upward(krho,idx-1); J=1i;
            Zr.Ze=Z.Ze.*(Zr_next.Ze+J.*Z.Ze.*t)./(Z.Ze+J.*Zr_next.Ze.*t); Zr.Zh=Z.Zh.*(Zr_next.Zh+J.*Z.Zh.*t)./(Z.Zh+J.*Zr_next.Zh.*t); Zr.Ye=1./Zr.Ze; Zr.Yh=1./Zr.Zh;
        end
        function Zl=ComputeInterfaceImpedance_Downward(obj, krho, idx)
            N=length(obj.lm.layers); if idx==N-1, Zl=obj.Z_bot; return; end
            % Fix: isPEC_bot
            if idx==N-2 && obj.lm.isPEC_bot, [~,~,k]=obj.GetLayerParams(idx+1); kz=obj.ComputeAxialWaveNumber(k,krho); t=tan(kz*obj.lm.GetHeight(idx+1)); Z=obj.ComputeLayerImpedance(krho,idx+1,kz); J=1i; Zl.Ze=Z.Ze.*J.*t; Zl.Zh=Z.Zh.*J.*t; Zl.Ye=Z.Ye./(J.*t); Zl.Yh=Z.Yh./(J.*t); return; end
            [~,~,k]=obj.GetLayerParams(idx+1); kz=obj.ComputeAxialWaveNumber(k,krho); t=tan(kz*obj.lm.GetHeight(idx+1)); Z=obj.ComputeLayerImpedance(krho,idx+1,kz); Zl_next=obj.ComputeInterfaceImpedance_Downward(krho,idx+1); J=1i;
            Zl.Ze=Z.Ze.*(Zl_next.Ze+J.*Z.Ze.*t)./(Z.Ze+J.*Zl_next.Ze.*t); Zl.Zh=Z.Zh.*(Zl_next.Zh+J.*Z.Zh.*t)./(Z.Zh+J.*Zl_next.Zh.*t); Zl.Ye=1./Zl.Ze; Zl.Yh=1./Zl.Zh;
        end
        % Reflection Coeffs (Fixes applied internally via ComputeHalfspaceImpedance isPEC checks)
        function Gr=ComputeReflectionCoefficient_Upward(obj, krho, idx)
            % Fix: isPEC_top
            if idx==-2
                Gr.Ge=0; 
                Gr.Gh=0; 
                return; 
            end
            if idx==-1 && obj.lm.isPEC_top 
                Gr.Ge=-1; 
                Gr.Gh=-1; 
                return; 
            end
            Z=obj.ComputeLayerImpedance(krho,idx+1,[]); 
            Zr=obj.ComputeInterfaceImpedance_Upward(krho,idx);
            Gr.Ge=(Zr.Ze-Z.Ze)./(Zr.Ze+Z.Ze); 
            Gr.Gh=(Zr.Zh-Z.Zh)./(Zr.Zh+Z.Zh);
        end
        function Gl=ComputeReflectionCoefficient_Downward(obj, krho, idx)
            N=length(obj.lm.layers); 
            if idx==N
                Gl.Ge=0; 
                Gl.Gh=0; 
                return; 
            end
            if idx==N-1 && obj.lm.isPEC_bot
                Gl.Ge=-1; 
                Gl.Gh=-1; 
                return; 
            end
            Z=obj.ComputeLayerImpedance(krho,idx,[]); 
            Zl=obj.ComputeInterfaceImpedance_Downward(krho,idx);
            Gl.Ge=(Zl.Ze-Z.Ze)./(Zl.Ze+Z.Ze); 
            Gl.Gh=(Zl.Zh-Z.Zh)./(Zl.Zh+Z.Zh);
        end
        function [MVe,MVh,MIe,MIh]=ComputeM_Upward(obj, krho)
            MVe=1; MVh=1; MIe=1; MIh=1; J=1i;
            for ii=obj.m:(obj.i-2)
                hk=obj.lm.GetHeight(ii+1); 
                [~,~,kk]=obj.GetLayerParams(ii+1); 
                kzk=obj.ComputeAxialWaveNumber(kk,krho); 
                jkzk=J*kzk;
                Gk=obj.ComputeReflectionCoefficient_Upward(krho,ii); 
                num_Ve=(1+Gk.Ge).*exp(-jkzk*hk); 
                den_Ve=(1+Gk.Ge.*exp(-2*jkzk*hk));                 
                MVe=MVe.*(num_Ve./den_Ve);

                num_Vh=(1+Gk.Gh).*exp(-jkzk*hk); 
                den_Vh=(1+Gk.Gh.*exp(-2*jkzk*hk)); 
                MVh=MVh.*(num_Vh./den_Vh);
                
                num_Ie = (1-Gk.Ge).*exp(-jkzk*hk);
                den_Ie = (1-Gk.Ge.*exp(-2*jkzk*hk));
                MIe=MIe.*(num_Ie./den_Ie); 

                num_Ih = (1-Gk.Gh).*exp(-jkzk*hk);
                den_Ih = (1-Gk.Gh.*exp(-2*jkzk*hk));
                MIh=MIh.*(num_Ih./den_Ih);
            end
            hk=obj.lm.GetHeight(obj.m); 
            kzk=obj.kzm; 
            jkzk=1i*kzk;
            MVe=MVe.*(1./(1+obj.Gm.Ge.*exp(-jkzk*2.0*hk))); 
            MVh=MVh.*(1./(1+obj.Gm.Gh.*exp(-jkzk*2.0*hk)));
            MIe=MIe.*(1./(1-obj.Gm.Ge.*exp(-jkzk*2.0*hk))); 
            MIh=MIh.*(1./(1-obj.Gm.Gh.*exp(-jkzk*2.0*hk)));
        end
        function [MVe,MVh,MIe,MIh]=ComputeM_Downward(obj, krho)
            MVe=1; MVh=1; MIe=1; MIh=1; J=1i;
            for ii=(obj.i+1):(obj.m-1)
                hk=obj.lm.GetHeight(ii); 
                [~,~,kk]=obj.GetLayerParams(ii); 
                kzk=obj.ComputeAxialWaveNumber(kk,krho); 
                jkzk=J*kzk;
                Gk=obj.ComputeReflectionCoefficient_Downward(krho,ii); 

                num_Ve=(1+Gk.Ge).*exp(-jkzk*hk); 
                den_Ve=(1+Gk.Ge.*exp(-2*jkzk*hk)); 
                MVe=MVe.*(num_Ve./den_Ve);

                num_Vh=(1+Gk.Gh).*exp(-jkzk*hk); 
                den_Vh=(1+Gk.Gh.*exp(-2*jkzk*hk)); 
                MVh=MVh.*(num_Vh./den_Vh);

                num_Ie = (1-Gk.Ge).*exp(-jkzk*hk);
                den_Ie = (1-Gk.Ge.*exp(-2*jkzk*hk));
                MIe=MIe.*(num_Ie./den_Ie); 

                num_Ih = (1-Gk.Gh).*exp(-jkzk*hk);
                den_Ih = (1-Gk.Gh.*exp(-2*jkzk*hk));
                MIh=MIh.*(num_Ih./den_Ih);
            end
            hk=obj.lm.GetHeight(obj.m); 
            kzk=obj.kzm; jkzk=1i*kzk;
            MVe=MVe.*(1./(1+obj.Gm.Ge.*exp(-jkzk*2.0*hk)));
            MVh=MVh.*(1./(1+obj.Gm.Gh.*exp(-jkzk*2.0*hk)));
            MIe=MIe.*(1./(1-obj.Gm.Ge.*exp(-jkzk*2.0*hk))); 
            MIh=MIh.*(1./(1-obj.Gm.Gh.*exp(-jkzk*2.0*hk)));
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

        function [eps,mu,k]=GetLayerParams(obj,idx)
            [eps,mu,k]=obj.lm.getLayerParams(idx); 
        end
        function ResetComputedFlags(obj)
            obj.exp_terms_computed=false; 
            obj.trf_terms_computed=false; 
            obj.src_terms_computed=false; 
            obj.obs_terms_computed=false; 
        end
        function LogVariables(obj)      % 我调试bug用的
            fid = fopen('debug_spectral.txt', 'a');
            if fid == -1, return; end
            
            % 格式化函数：复数转 (re,im) 字符串
            fmt = @(x) sprintf('(%0.5e,%0.5e)', real(x), imag(x));
            
            % 提取需要打印的变量 (假设 ComputeKmi 中已计算)
            str = sprintf('TIe: %s, TIh: %s, GVe_ii: %s, GVh_ii: %s, GIe_ii: %s, GIh_ii: %s, ki: %s, krhosq: %s, kzm: %s, km: %s, K[3]: %s, K[2]: %s\n', ...
                fmt(obj.TIe), fmt(obj.TIh), ...
                fmt(obj.GVe), fmt(obj.GVh), ...
                fmt(obj.GIe), fmt(obj.GIh), ...
                fmt(obj.ki), fmt(obj.krhosq), ...
                fmt(obj.kzm), fmt(obj.km), ...
                fmt(obj.K(4)), fmt(obj.K(3)));
            
            fprintf(fid, '%s', str);
            fclose(fid);
        end
    end
end