%% 最终方案：用 k_max 作为 DCIM 路径参数
clear classes; %#ok<CLCLS>
clear; clc;
addpath(pwd);

freq = 30e9; unit = 1e-3;
c0 = 1/sqrt(8.854187817e-12 * 4*pi*1e-7);
lambda0 = c0/freq;

lm = LayerManager();
lm.AddLayer(1.1*unit, 1.8*unit, 2.1, 1, 0);
lm.AddLayer(0.8*unit, 1.1*unit, 12.5, 1, 0);
lm.AddLayer(0.3*unit, 0.8*unit, 9.8, 1, 0);
lm.AddLayer(0.0*unit, 0.3*unit, 8.6, 1, 0);
lm.SetHalfspaces(1.0, 1.0, 0.0, 1.0, 1.0, 0.0, false, true);
lm.ProcessLayers(freq);

z_src = 0.4*unit; z_obs = 1.4*unit;
i_src = lm.FindLayer(z_src); m_obs = lm.FindLayer(z_obs);
rho_test = 1e-4 * lambda0;

% 基准
smgf2 = SpectralMGF();
smgf2.Initialize(lm, freq, [1,1,1,1,1], true, false);
smgf2.SetLayers(i_src, m_obs);
smgf2.SetSourcePoint(z_src); smgf2.SetObservationPoint(z_obs);
a_sw = 1.2*abs(lm.k_max);
qmgf = QuasistaticMGF2();
qmgf.Initialize(lm, freq, [1,1,1,1,1], false);
qmgf.SetLayers(i_src, m_obs);
K_q = qmgf.ComputeQMGF_Spatial(z_obs, z_src, rho_test);
orders = [0,1,1,0,0];
G_ref = zeros(1,5);
for cc=1:5
    G_int = SommerfeldIntegrator2.IntegrateSpectralNearField(smgf2, rho_test, cc, orders(cc), a_sw, false) ...
          + SommerfeldIntegrator2.IntegrateSpectralFarField(smgf2, rho_test, cc, orders(cc), a_sw, false);
    G_ref(cc) = G_int + K_q(cc);
end

epsr_list = [lm.layers.epsr]; epsr_max = max(epsr_list);

% 扫描 k 和 T02 对所有分量的综合误差
fprintf('=== 扫描 k 和 T02（综合所有分量） ===\n');
best_total_err = inf;
for k_mult = [1.0, 1.2, 1.5, 1.8, 2.0, 2.44]
    kk = lm.k_min * k_mult;
    for T02_f = [5, 8, 12, 20, 30, 40]
        T02 = T02_f * sqrt(epsr_max);
        
        errs = zeros(1,5);
        for comp = 1:5
            errs(comp) = run_dcim(lm, freq, i_src, m_obs, z_src, z_obs, rho_test, kk, T02, comp, G_ref(comp));
        end
        
        % 加权误差：主对角分量(1,4,5)权重高
        total_err = (errs(1) + errs(4) + errs(5)) / 3;
        
        if total_err < best_total_err
            best_total_err = total_err;
            best_k = kk;
            best_T02 = T02;
            best_errs = errs;
        end
        
        fprintf('k=%7.1f(%.1fx) T02=%5.1f(%2dx): Gxx=%5.1f%% Gxz=%5.1f%% Gzx=%5.1f%% Gzz=%5.1f%% Gphi=%5.1f%% avg_diag=%5.1f%%\n', ...
            kk, k_mult, T02, T02_f, errs(1), errs(2), errs(3), errs(4), errs(5), total_err);
    end
end

fprintf('\n最优: k=%.1f (%.2f*k_min), T02=%.1f (%d*sqrt(epsr_max))\n', ...
    best_k, best_k/lm.k_min, best_T02, round(best_T02/sqrt(epsr_max)));
fprintf('误差: Gxx=%.2f%% Gxz=%.2f%% Gzx=%.2f%% Gzz=%.2f%% Gphi=%.2f%%\n', ...
    best_errs(1), best_errs(2), best_errs(3), best_errs(4), best_errs(5));

function err = run_dcim(lm, freq, i_src, m_obs, z_src, z_obs, rho_test, kk, T02, comp, G_ref_val)
    T01 = 200.0; N1 = 100; N2 = 100;
    
    smgf = SpectralMGF();
    smgf.Initialize(lm, freq, [1,1,1,1,1], false, false);
    smgf.SetLayers(i_src, m_obs);
    smgf.SetSourcePoint(z_src); smgf.SetObservationPoint(z_obs);
    
    t1 = linspace(0, T01, N1);
    kz1 = -1j*kk*(T02+t1);
    krho1 = sqrt(kk^2-kz1.^2); krho1 = abs(real(krho1))+1j*abs(imag(krho1));
    t2 = linspace(1e-2, T02, N2);
    kz2 = kk*(-1j*t2+(1-t2/T02));
    krho2 = sqrt(kk^2-kz2.^2); krho2 = abs(real(krho2))+1j*abs(imag(krho2));
    
    K1 = zeros(1,N1); K2 = zeros(1,N2);
    for ii=1:N1
        smgf.SetRadialWaveNumber(krho1(ii)); smgf.ComputeSpectralMGF();
        K1(ii) = smgf.K(comp)*kz1(ii);
    end
    for ii=1:N2
        smgf.SetRadialWaveNumber(krho2(ii)); smgf.ComputeSpectralMGF();
        K2(ii) = smgf.K(comp)*kz2(ii);
    end
    
    [at1, a1, ~] = GPOF(K1(:), t1(:), 50, 1e-4, 1e-16);
    a1k = a1.*exp(-T02*at1); al1k = at1./(1j*kk);
    
    K2m = K2;
    for ii=1:N2
        for jj=1:length(a1k)
            K2m(ii) = K2m(ii)-a1k(jj)*exp(-al1k(jj)*kz2(ii));
        end
    end
    
    [at2, a2, ~] = GPOF(K2m(:), t2(:), 50, 1e-4, 1e-16);
    al2k = at2.*T02./((1+1j*T02)*kk);
    a2k = a2.*exp(kk*al2k);
    
    aa = [a1k(:);a2k(:)]; aal = [al1k(:);al2k(:)];
    
    G = 0; rho_sq = rho_test^2;
    for ii=1:length(aa)
        Rc = sqrt(rho_sq-aal(ii)^2);
        if comp==1||comp==4||comp==5
            G = G + aa(ii)*exp(-1j*kk*Rc)/Rc;
        else
            G = G + rho_test*aa(ii)*(1+1j*kk*Rc)*exp(-1j*kk*Rc)/(Rc^3);
        end
    end
    G = G*1j/(2*pi);
    
    err = abs(G-G_ref_val)/abs(G_ref_val)*100;
end
