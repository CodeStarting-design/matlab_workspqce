%% 诊断准静态减法 + 尝试三层DCIM替代方案
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
comp = 1;

fprintf('k = %.4f, T02 = %.4f\n', k, T02);

%% 检查准静态谱域减法
smgf_full = SpectralMGF();
smgf_full.Initialize(lm, freq, [1,1,1,1,1], false, false);
smgf_full.SetLayers(i_src, m_obs);
smgf_full.SetSourcePoint(z_src);
smgf_full.SetObservationPoint(z_obs);

smgf_qs = SpectralMGF();
smgf_qs.Initialize(lm, freq, [1,1,1,1,1], true, false);
smgf_qs.SetLayers(i_src, m_obs);
smgf_qs.SetSourcePoint(z_src);
smgf_qs.SetObservationPoint(z_obs);

fprintf('\n=== 谱域函数在路径2端点处 ===\n');
N2 = 100;
t2 = linspace(1e-2, T02, N2);
kz2 = k*(-1j*t2+(1-t2/T02));
krho2 = sqrt(k^2-kz2.^2);
krho2 = abs(real(krho2))+1j*abs(imag(krho2));

K_full = zeros(1,N2);
K_qs = zeros(1,N2);
for ii = 1:N2
    smgf_full.SetRadialWaveNumber(krho2(ii));
    smgf_full.ComputeSpectralMGF();
    K_full(ii) = smgf_full.K(comp) * kz2(ii);
    
    smgf_qs.SetRadialWaveNumber(krho2(ii));
    smgf_qs.ComputeSpectralMGF();
    K_qs(ii) = smgf_qs.K(comp) * kz2(ii);
end

fprintf('Full:  max=%.4e, norm=%.4e\n', max(abs(K_full)), norm(K_full));
fprintf('QS extracted: max=%.4e, norm=%.4e\n', max(abs(K_qs)), norm(K_qs));
fprintf('Ratio(norm): %.2f\n', norm(K_qs)/norm(K_full));

%% 三层DCIM方案 (参考strata DCIM.cpp GenerateSamplePoints_ThreeLevel)
fprintf('\n======== 三层 DCIM ========\n');

T03 = 2.0*sqrt(epsr_max);  % strata ThreeLevel
T02_3 = 200.0;
T01_3 = 2000.0;
N1 = 100; N2 = 100; N3 = 100;

smgf = SpectralMGF();
smgf.Initialize(lm, freq, [1,1,1,1,1], false, false);
smgf.SetLayers(i_src, m_obs);
smgf.SetSourcePoint(z_src);
smgf.SetObservationPoint(z_obs);

% Path 1: kz = -j*k*(T03+T02+t), t in [0, T01]
t1 = linspace(0, T01_3, N1);
kz1 = -1j*k*(T03+T02_3+t1);
krho1 = sqrt(k^2-kz1.^2);
krho1 = abs(real(krho1))+1j*abs(imag(krho1));

% Path 2: kz = -j*k*(T03+t), t in [0, T02]
t2 = linspace(0, T02_3, N2);
kz2 = -1j*k*(T03+t2);
krho2 = sqrt(k^2-kz2.^2);
krho2 = abs(real(krho2))+1j*abs(imag(krho2));

% Path 3: kz = k*(-j*t + (1 - t/T03)), t in [1e-8, T03]
t3 = linspace(1e-8, T03, N3);
kz3 = k*(-1j*t3+(1-t3/T03));
krho3 = sqrt(k^2-kz3.^2);
krho3 = abs(real(krho3))+1j*abs(imag(krho3));

% 采样 path 1
K1 = zeros(1,N1);
for ii=1:N1
    smgf.SetRadialWaveNumber(krho1(ii));
    smgf.ComputeSpectralMGF();
    K1(ii) = smgf.K(comp)*kz1(ii);
end

% 采样 path 2
K2 = zeros(1,N2);
for ii=1:N2
    smgf.SetRadialWaveNumber(krho2(ii));
    smgf.ComputeSpectralMGF();
    K2(ii) = smgf.K(comp)*kz2(ii);
end

% 采样 path 3
K3 = zeros(1,N3);
for ii=1:N3
    smgf.SetRadialWaveNumber(krho3(ii));
    smgf.ComputeSpectralMGF();
    K3(ii) = smgf.K(comp)*kz3(ii);
end

fprintf('Path1: max=%.4e\n', max(abs(K1)));
fprintf('Path2: max=%.4e\n', max(abs(K2)));
fprintf('Path3: max=%.4e\n', max(abs(K3)));

tol_svd = 1e-4; tol_eig = 1e-16;

% Level 1 GPOF
L1 = floor(N1/2);
[at1, a1, res1] = GPOF(K1(:), t1(:), L1, tol_svd, tol_eig);
fprintf('\nL1: %d exp, res=%.4e\n', length(at1), res1);

% t->kz: a *= exp(-(T03+T02)*alpha), alpha /= (j*k)
a1k = a1.*exp(-(T03+T02_3)*at1);
al1k = at1./(1j*k);

% Subtract L1 from path2 and path3
K2m = K2;
for ii=1:N2
    for jj=1:length(a1k)
        K2m(ii) = K2m(ii)-a1k(jj)*exp(-al1k(jj)*kz2(ii));
    end
end
K3m = K3;
for ii=1:N3
    for jj=1:length(a1k)
        K3m(ii) = K3m(ii)-a1k(jj)*exp(-al1k(jj)*kz3(ii));
    end
end

% Level 2 GPOF
L2 = floor(N2/2);
[at2, a2, res2] = GPOF(K2m(:), t2(:), L2, tol_svd, tol_eig);
fprintf('L2: %d exp, res=%.4e\n', length(at2), res2);

% t->kz: a *= exp(-T03*alpha), alpha /= (j*k)
a2k = a2.*exp(-T03*at2);
al2k = at2./(1j*k);

% Subtract L2 from path3
for ii=1:N3
    for jj=1:length(a2k)
        K3m(ii) = K3m(ii)-a2k(jj)*exp(-al2k(jj)*kz3(ii));
    end
end

% Level 3 GPOF
L3 = floor(N3/2);
[at3, a3, res3] = GPOF(K3m(:), t3(:), L3, tol_svd, tol_eig);
fprintf('L3: %d exp, res=%.4e\n', length(at3), res3);

% t->kz: alpha = alpha*T03/((1+j*T03)*k), a = a*exp(k*alpha)
al3k = at3.*T03./((1+1j*T03)*k);
a3k = a3.*exp(k*al3k);

% 合并
aa = [a1k(:);a2k(:);a3k(:)];
aal = [al1k(:);al2k(:);al3k(:)];
fprintf('总镜像: %d\n', length(aa));

% 空间域
G_dcim = 0;
for ii=1:length(aa)
    Rc = sqrt(rho_test^2-aal(ii)^2);
    G_dcim = G_dcim + aa(ii)*exp(-1j*k*Rc)/Rc;
end
G_dcim = G_dcim*1j/(2*pi);

% 基准
smgf2 = SpectralMGF();
smgf2.Initialize(lm, freq, [1,1,1,1,1], true, false);
smgf2.SetLayers(i_src, m_obs);
smgf2.SetSourcePoint(z_src); smgf2.SetObservationPoint(z_obs);
a_sw = 1.2*abs(lm.k_max);
G_int = SommerfeldIntegrator2.IntegrateSpectralNearField(smgf2, rho_test, comp, 0, a_sw, false) ...
      + SommerfeldIntegrator2.IntegrateSpectralFarField(smgf2, rho_test, comp, 0, a_sw, false);
qmgf2 = QuasistaticMGF2();
qmgf2.Initialize(lm, freq, [1,1,1,1,1], false);
qmgf2.SetLayers(i_src, m_obs);
K_q2 = qmgf2.ComputeQMGF_Spatial(z_obs, z_src, rho_test);
G_ref = G_int + K_q2(comp);

fprintf('\n三层DCIM:\n');
fprintf('G_dcim = %.6e %+.6ei\n', real(G_dcim), imag(G_dcim));
fprintf('G_ref  = %.6e %+.6ei\n', real(G_ref), imag(G_ref));
fprintf('误差 = %.2f%%\n', abs(G_dcim-G_ref)/abs(G_ref)*100);

% 实轴验证
fprintf('\n实轴验证:\n');
krho_tests = [10, 500, 800, 1000, 1500, 2000];
for kk = 1:length(krho_tests)
    kr = krho_tests(kk)+0.001j;
    kz_t = sqrt(k^2-kr^2);
    if imag(kz_t)>0, kz_t = complex(real(kz_t),-abs(imag(kz_t))); end
    if real(kz_t)<0, kz_t = complex(abs(real(kz_t)),imag(kz_t)); end
    smgf.SetRadialWaveNumber(kr);
    smgf.ComputeSpectralMGF();
    K_dir = smgf.K(comp)*kz_t;
    K_dcim = 0;
    for jj=1:length(aa)
        K_dcim = K_dcim+aa(jj)*exp(-aal(jj)*kz_t);
    end
    fprintf('  krho=%5.0f: err=%.1f%%\n', real(kr), abs(K_dcim-K_dir)/abs(K_dir)*100);
end
