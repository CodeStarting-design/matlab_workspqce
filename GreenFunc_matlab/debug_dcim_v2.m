%% 诊断准静态一致性 + 测试 extract_quasistatic=true 的 DCIM
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
comp = 1;
fprintf('k_min = %.4f\n', k);

%% 第1步：验证准静态一致性
% 检查 SpectralMGF(qs=true) + QuasistaticSpatial == SpectralMGF(qs=false)
% 即谱域中：K_full = K_qs_extracted + K_quasi_spectral
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

qmgf_spec = QuasistaticMGF2();
qmgf_spec.Initialize(lm, freq, [1,1,1,1,1], false);
qmgf_spec.SetLayers(i_src, m_obs);

fprintf('\n=== 准静态谱域一致性验证 ===\n');
krho_checks = [10, 100, 500, 900, 1000, 1500, 2000];
for kk = 1:length(krho_checks)
    kr = krho_checks(kk) + 0.001j;
    
    smgf_full.SetRadialWaveNumber(kr);
    smgf_full.ComputeSpectralMGF();
    K_full = smgf_full.K(comp);
    
    smgf_qs.SetRadialWaveNumber(kr);
    smgf_qs.ComputeSpectralMGF();
    K_extracted = smgf_qs.K(comp);
    
    K_quasi_spec = qmgf_spec.ComputeQMGF_Spectral(z_obs, z_src, kr);
    K_quasi_val = K_quasi_spec(comp);
    
    K_recon = K_extracted + K_quasi_val;
    err_recon = abs(K_recon - K_full)/abs(K_full)*100;
    
    fprintf('krho=%5.0f: |K_full|=%.3e, |K_ext|=%.3e, |K_qs|=%.3e, recon_err=%.2f%%\n', ...
        real(kr), abs(K_full), abs(K_extracted), abs(K_quasi_val), err_recon);
end

%% 第2步：用 extract_quasistatic=true 做两层 DCIM
fprintf('\n=== 两层 DCIM with extract_quasistatic=true ===\n');
epsr_list = [lm.layers.epsr]; epsr_max = max(epsr_list);
T02 = 5.0*sqrt(epsr_max); T01 = 200.0;
N1 = 100; N2 = 100;

t1 = linspace(0, T01, N1);
kz1 = -1j*k*(T02+t1);
krho1 = sqrt(k^2-kz1.^2);
krho1 = abs(real(krho1))+1j*abs(imag(krho1));

t2 = linspace(1e-2, T02, N2);
kz2 = k*(-1j*t2+(1-t2/T02));
krho2 = sqrt(k^2-kz2.^2);
krho2 = abs(real(krho2))+1j*abs(imag(krho2));

% 采样（使用已提取准静态的 smgf_qs）
K1 = zeros(1,N1);
for ii=1:N1
    smgf_qs.SetRadialWaveNumber(krho1(ii));
    smgf_qs.ComputeSpectralMGF();
    K1(ii) = smgf_qs.K(comp)*kz1(ii);
end

K2 = zeros(1,N2);
for ii=1:N2
    smgf_qs.SetRadialWaveNumber(krho2(ii));
    smgf_qs.ComputeSpectralMGF();
    K2(ii) = smgf_qs.K(comp)*kz2(ii);
end

fprintf('Path1: max=%.4e, norm=%.4e\n', max(abs(K1)), norm(K1));
fprintf('Path2: max=%.4e, norm=%.4e\n', max(abs(K2)), norm(K2));

% Level 1 GPOF
[at1, a1, res1] = GPOF(K1(:), t1(:), 50, 1e-4, 1e-16);
fprintf('L1: %d exp, res=%.4e\n', length(at1), res1);
a1k = a1.*exp(-T02*at1);
al1k = at1./(1j*k);

% Subtract L1 from path2
K2m = K2;
for ii=1:N2
    for jj=1:length(a1k)
        K2m(ii) = K2m(ii)-a1k(jj)*exp(-al1k(jj)*kz2(ii));
    end
end
fprintf('After L1 subtract: max|K2m|=%.4e (was %.4e)\n', max(abs(K2m)), max(abs(K2)));

% Level 2 GPOF
[at2, a2, res2] = GPOF(K2m(:), t2(:), 50, 1e-4, 1e-16);
fprintf('L2: %d exp, res=%.4e\n', length(at2), res2);
al2k = at2.*T02./((1+1j*T02)*k);
a2k = a2.*exp(k*al2k);

aa = [a1k(:);a2k(:)];
aal = [al1k(:);al2k(:)];
fprintf('总镜像: %d\n', length(aa));

% 空间域
G_dcim_core = 0;
for ii=1:length(aa)
    Rc = sqrt(rho_test^2-aal(ii)^2);
    G_dcim_core = G_dcim_core + aa(ii)*exp(-1j*k*Rc)/Rc;
end
G_dcim_core = G_dcim_core*1j/(2*pi);

% 加回空间域准静态
qmgf_sp = QuasistaticMGF2();
qmgf_sp.Initialize(lm, freq, [1,1,1,1,1], false);
qmgf_sp.SetLayers(i_src, m_obs);
K_quasi_spatial = qmgf_sp.ComputeQMGF_Spatial(z_obs, z_src, rho_test);
G_dcim = G_dcim_core + K_quasi_spatial(comp);

% 基准
smgf_ref = SpectralMGF();
smgf_ref.Initialize(lm, freq, [1,1,1,1,1], true, false);
smgf_ref.SetLayers(i_src, m_obs);
smgf_ref.SetSourcePoint(z_src);
smgf_ref.SetObservationPoint(z_obs);
a_sw = 1.2*abs(lm.k_max);
G_int = SommerfeldIntegrator2.IntegrateSpectralNearField(smgf_ref, rho_test, comp, 0, a_sw, false) ...
      + SommerfeldIntegrator2.IntegrateSpectralFarField(smgf_ref, rho_test, comp, 0, a_sw, false);
G_ref = G_int + K_quasi_spatial(comp);

fprintf('\nG_dcim_core = %.6e %+.6ei\n', real(G_dcim_core), imag(G_dcim_core));
fprintf('G_int_core  = %.6e %+.6ei\n', real(G_int), imag(G_int));
fprintf('K_quasi_sp  = %.6e %+.6ei\n', real(K_quasi_spatial(comp)), imag(K_quasi_spatial(comp)));
fprintf('G_dcim      = %.6e %+.6ei\n', real(G_dcim), imag(G_dcim));
fprintf('G_ref       = %.6e %+.6ei\n', real(G_ref), imag(G_ref));
fprintf('误差 = %.2f%%\n', abs(G_dcim - G_ref)/abs(G_ref)*100);

% 实轴验证（提取后）
fprintf('\n实轴验证（提取后 DCIM 核心）:\n');
krho_tests = [10, 500, 800, 1000, 1500, 2000];
for kk = 1:length(krho_tests)
    kr = krho_tests(kk)+0.001j;
    kz_t = sqrt(k^2-kr^2);
    if imag(kz_t)>0, kz_t = complex(real(kz_t),-abs(imag(kz_t))); end
    if real(kz_t)<0, kz_t = complex(abs(real(kz_t)),imag(kz_t)); end
    smgf_qs.SetRadialWaveNumber(kr);
    smgf_qs.ComputeSpectralMGF();
    K_dir = smgf_qs.K(comp)*kz_t;
    K_dcim_v = 0;
    for jj=1:length(aa)
        K_dcim_v = K_dcim_v+aa(jj)*exp(-aal(jj)*kz_t);
    end
    fprintf('  krho=%5.0f: |K_dir|=%.3e, |K_dcim|=%.3e, err=%.1f%%\n', ...
        real(kr), abs(K_dir), abs(K_dcim_v), abs(K_dcim_v-K_dir)/max(abs(K_dir),1e-30)*100);
end
