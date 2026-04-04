%% 扫描 T02 参数优化
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
k = lm.k_min;

% 基准
smgf2 = SpectralMGF();
smgf2.Initialize(lm, freq, [1,1,1,1,1], true, false);
smgf2.SetLayers(i_src, m_obs);
smgf2.SetSourcePoint(z_src); smgf2.SetObservationPoint(z_obs);
a_sw = 1.2*abs(lm.k_max);
G_int = SommerfeldIntegrator2.IntegrateSpectralNearField(smgf2, rho_test, 1, 0, a_sw, false) ...
      + SommerfeldIntegrator2.IntegrateSpectralFarField(smgf2, rho_test, 1, 0, a_sw, false);
qmgf = QuasistaticMGF2();
qmgf.Initialize(lm, freq, [1,1,1,1,1], false);
qmgf.SetLayers(i_src, m_obs);
K_q = qmgf.ComputeQMGF_Spatial(z_obs, z_src, rho_test);
G_ref = G_int + K_q(1);
fprintf('G_ref = abs %.4f\n\n', abs(G_ref));

% 用 k = k_max 代替 k_min 试试（strata的epsr_max计算可能影响k选择）
fprintf('k_min=%.2f, k_max=%.2f\n\n', lm.k_min, abs(lm.k_max));

for k_choice = [lm.k_min, abs(lm.k_max), (lm.k_min+abs(lm.k_max))/2]
    for T02_factor = [3, 5, 8, 12]
        epsr_list = [lm.layers.epsr]; epsr_max = max(epsr_list);
        T02 = T02_factor * sqrt(epsr_max);
        T01 = 200.0;
        N1 = 100; N2 = 100;
        kk = k_choice;
        
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
            K1(ii) = smgf.K(1)*kz1(ii);
        end
        for ii=1:N2
            smgf.SetRadialWaveNumber(krho2(ii)); smgf.ComputeSpectralMGF();
            K2(ii) = smgf.K(1)*kz2(ii);
        end
        
        [at1, a1, ~] = GPOF(K1(:), t1(:), 50, 1e-4, 1e-16);
        a1k = a1.*exp(-T02*at1); al1k = at1./(1j*kk);
        
        K2m = K2;
        for ii=1:N2
            for jj=1:length(a1k)
                K2m(ii) = K2m(ii)-a1k(jj)*exp(-al1k(jj)*kz2(ii));
            end
        end
        
        [at2, a2, res2] = GPOF(K2m(:), t2(:), 50, 1e-4, 1e-16);
        al2k = at2.*T02./((1+1j*T02)*kk);
        a2k = a2.*exp(kk*al2k);
        
        aa = [a1k(:);a2k(:)]; aal = [al1k(:);al2k(:)];
        G = compute_G(aa, aal, kk, rho_test);
        err = abs(G-G_ref)/abs(G_ref)*100;
        fprintf('k=%7.1f T02=%5.1f | L2:%2d res=%.3e max|a|=%.2e | err=%6.2f%%\n', ...
            kk, T02, length(at2), res2, max(abs(a2k)), err);
    end
    fprintf('\n');
end

function G = compute_G(aa, aal, k, rho)
    G = 0;
    for ii=1:length(aa)
        Rc = sqrt(rho^2-aal(ii)^2);
        G = G + aa(ii)*exp(-1j*k*Rc)/Rc;
    end
    G = G*1j/(2*pi);
end
