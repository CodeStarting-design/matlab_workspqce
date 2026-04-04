%% 测试镜像过滤策略
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
epsr_list = [lm.layers.epsr]; epsr_max = max(epsr_list);
T02 = 5.0*sqrt(epsr_max); T01 = 200.0;
N1 = 100; N2 = 100; comp = 1;

smgf = SpectralMGF();
smgf.Initialize(lm, freq, [1,1,1,1,1], false, false);
smgf.SetLayers(i_src, m_obs);
smgf.SetSourcePoint(z_src); smgf.SetObservationPoint(z_obs);

t1 = linspace(0, T01, N1);
kz1 = -1j*k*(T02+t1);
krho1 = sqrt(k^2-kz1.^2); krho1 = abs(real(krho1))+1j*abs(imag(krho1));
t2 = linspace(1e-2, T02, N2);
kz2 = k*(-1j*t2+(1-t2/T02));
krho2 = sqrt(k^2-kz2.^2); krho2 = abs(real(krho2))+1j*abs(imag(krho2));

K1 = zeros(1,N1);
for ii=1:N1
    smgf.SetRadialWaveNumber(krho1(ii)); smgf.ComputeSpectralMGF();
    K1(ii) = smgf.K(comp)*kz1(ii);
end
K2 = zeros(1,N2);
for ii=1:N2
    smgf.SetRadialWaveNumber(krho2(ii)); smgf.ComputeSpectralMGF();
    K2(ii) = smgf.K(comp)*kz2(ii);
end

% 基准
smgf2 = SpectralMGF();
smgf2.Initialize(lm, freq, [1,1,1,1,1], true, false);
smgf2.SetLayers(i_src, m_obs);
smgf2.SetSourcePoint(z_src); smgf2.SetObservationPoint(z_obs);
a_sw = 1.2*abs(lm.k_max);
G_int = SommerfeldIntegrator2.IntegrateSpectralNearField(smgf2, rho_test, comp, 0, a_sw, false) ...
      + SommerfeldIntegrator2.IntegrateSpectralFarField(smgf2, rho_test, comp, 0, a_sw, false);
qmgf = QuasistaticMGF2();
qmgf.Initialize(lm, freq, [1,1,1,1,1], false);
qmgf.SetLayers(i_src, m_obs);
K_q = qmgf.ComputeQMGF_Spatial(z_obs, z_src, rho_test);
G_ref = G_int + K_q(comp);
fprintf('G_ref = abs %.4f\n\n', abs(G_ref));

% 用 tol_svd=1e-4（最优的）
tol_svd = 1e-4;
[at1, a1, ~] = GPOF(K1(:), t1(:), 50, tol_svd, 1e-16);
a1k = a1.*exp(-T02*at1); al1k = at1./(1j*k);

K2m = K2;
for ii=1:N2
    for jj=1:length(a1k)
        K2m(ii) = K2m(ii)-a1k(jj)*exp(-al1k(jj)*kz2(ii));
    end
end

[at2, a2, res2] = GPOF(K2m(:), t2(:), 50, tol_svd, 1e-16);
al2k = at2.*T02./((1+1j*T02)*k);
a2k = a2.*exp(k*al2k);

aa_all = [a1k(:);a2k(:)]; aal_all = [al1k(:);al2k(:)];

% 策略1: 不过滤（原始）
G1 = compute_G(aa_all, aal_all, k, rho_test);
fprintf('策略1 (原始 %d镜像): abs=%.4f, err=%.2f%%\n', length(aa_all), abs(G1), abs(G1-G_ref)/abs(G_ref)*100);

% 策略2: 过滤掉 |a_kz| > threshold 的镜像
for thresh = [50, 100, 500, 1000]
    mask = abs(aa_all) <= thresh;
    G2 = compute_G(aa_all(mask), aal_all(mask), k, rho_test);
    fprintf('策略2 (|a|<=%4d, %d镜像): abs=%.4f, err=%.2f%%\n', thresh, sum(mask), abs(G2), abs(G2-G_ref)/abs(G_ref)*100);
end

% 策略3: 用更多采样点 N2=200, 但限制 L2
fprintf('\n--- 更多采样点 ---\n');
for N2_test = [200, 300, 500]
    t2b = linspace(1e-2, T02, N2_test);
    kz2b = k*(-1j*t2b+(1-t2b/T02));
    krho2b = sqrt(k^2-kz2b.^2); krho2b = abs(real(krho2b))+1j*abs(imag(krho2b));
    
    K2b = zeros(1,N2_test);
    for ii=1:N2_test
        smgf.SetRadialWaveNumber(krho2b(ii)); smgf.ComputeSpectralMGF();
        K2b(ii) = smgf.K(comp)*kz2b(ii);
    end
    
    K2bm = K2b;
    for ii=1:N2_test
        for jj=1:length(a1k)
            K2bm(ii) = K2bm(ii)-a1k(jj)*exp(-al1k(jj)*kz2b(ii));
        end
    end
    
    % 用不同的L值
    for L2_test = [50, 80, 100]
        if L2_test >= N2_test, continue; end
        [at2b, a2b, res2b] = GPOF(K2bm(:), t2b(:), L2_test, tol_svd, 1e-16);
        al2kb = at2b.*T02./((1+1j*T02)*k);
        a2kb = a2b.*exp(k*al2kb);
        
        aab = [a1k(:);a2kb(:)]; aalb = [al1k(:);al2kb(:)];
        Gb = compute_G(aab, aalb, k, rho_test);
        fprintf('N2=%3d L2=%3d: %2d exp, res=%.3e, abs=%.4f, err=%.2f%%\n', ...
            N2_test, L2_test, length(at2b), res2b, abs(Gb), abs(Gb-G_ref)/abs(G_ref)*100);
    end
end

function G = compute_G(aa, aal, k, rho)
    G = 0;
    for ii=1:length(aa)
        Rc = sqrt(rho^2-aal(ii)^2);
        G = G + aa(ii)*exp(-1j*k*Rc)/Rc;
    end
    G = G*1j/(2*pi);
end
