%% 精确诊断 DCIM 22% 误差来源
clear classes; %#ok<CLCLS>
clear; clc;
addpath(pwd);

freq = 30e9;
unit = 1e-3;
c0 = 1/sqrt(8.854187817e-12 * 4*pi*1e-7);
lambda0 = c0/freq;

lm = LayerManager();
lm.AddLayer(1.1*unit, 1.8*unit, 2.1, 1, 0);
lm.AddLayer(0.8*unit, 1.1*unit, 12.5, 1, 0);
lm.AddLayer(0.3*unit, 0.8*unit, 9.8, 1, 0);
lm.AddLayer(0.0*unit, 0.3*unit, 8.6, 1, 0);
lm.SetHalfspaces(1.0, 1.0, 0.0, 1.0, 1.0, 0.0, false, true);
lm.ProcessLayers(freq);

z_src = 0.4 * unit;
z_obs = 1.4 * unit;
i_src = lm.FindLayer(z_src);
m_obs = lm.FindLayer(z_obs);
rho_test = 1e-4 * lambda0;

fprintf('rho_test = %.6e m = %.4f * lambda0\n', rho_test, rho_test/lambda0);
fprintf('k_min = %.4f\n', lm.k_min);

%% DCIM 手动 (extract_quasistatic=false)
k = lm.k_min;
epsr_list = [lm.layers.epsr];
epsr_max = max(epsr_list);
T02 = 5.0 * sqrt(epsr_max);
T01 = 200.0;
N1 = 100; N2 = 100;

smgf = SpectralMGF();
smgf.Initialize(lm, freq, [1,1,1,1,1], false, false);
smgf.SetLayers(i_src, m_obs);
smgf.SetSourcePoint(z_src);
smgf.SetObservationPoint(z_obs);

% 路径
t1 = linspace(0, T01, N1);
kz1 = -1j * k * (T02 + t1);
krho1 = sqrt(k^2 - kz1.^2);
krho1 = abs(real(krho1)) + 1j * abs(imag(krho1));

t2 = linspace(1e-2, T02, N2);
kz2 = k * (-1j * t2 + (1 - t2/T02));
krho2 = sqrt(k^2 - kz2.^2);
krho2 = abs(real(krho2)) + 1j * abs(imag(krho2));

comp = 1; % Gxx

% 路径1采样
K1 = zeros(1,N1);
for ii = 1:N1
    smgf.SetRadialWaveNumber(krho1(ii));
    smgf.ComputeSpectralMGF();
    K1(ii) = smgf.K(comp) * kz1(ii);
end

% GPOF L1
[alpha_t1, a1, res1] = GPOF(K1(:), t1(:), 50, 1e-4, 1e-16);
fprintf('\nL1: %d exp, res=%.4e\n', length(alpha_t1), res1);
a1_kz = a1 .* exp(-T02 * alpha_t1);
alpha1_kz = alpha_t1 ./ (1j * k);

% 路径2采样
K2 = zeros(1,N2);
for ii = 1:N2
    smgf.SetRadialWaveNumber(krho2(ii));
    smgf.ComputeSpectralMGF();
    K2(ii) = smgf.K(comp) * kz2(ii);
end

% 减去L1
K2_mod = K2;
for ii = 1:N2
    for jj = 1:length(a1_kz)
        K2_mod(ii) = K2_mod(ii) - a1_kz(jj) * exp(-alpha1_kz(jj) * kz2(ii));
    end
end
fprintf('L2 before GPOF: max|K2|=%.4e, max|K2_mod|=%.4e\n', max(abs(K2)), max(abs(K2_mod)));

% GPOF L2
[alpha_t2, a2, res2] = GPOF(K2_mod(:), t2(:), 50, 1e-4, 1e-16);
fprintf('L2: %d exp, res=%.4e\n', length(alpha_t2), res2);
alpha2_kz = alpha_t2 .* T02 ./ ((1 + 1j*T02) * k);
a2_kz = a2 .* exp(k * alpha2_kz);

all_a = [a1_kz(:); a2_kz(:)];
all_alpha = [alpha1_kz(:); alpha2_kz(:)];
fprintf('总复镜像数: %d\n', length(all_a));

%% 验证：在路径2上重建
K2_recon = zeros(1,N2);
for ii = 1:N2
    for jj = 1:length(all_a)
        K2_recon(ii) = K2_recon(ii) + all_a(jj) * exp(-all_alpha(jj) * kz2(ii));
    end
end
fprintf('路径2重建误差: %.4e\n', norm(K2-K2_recon)/norm(K2));

%% 验证：在实轴上重建
fprintf('\n=== 实轴验证 ===\n');
krho_tests = [10, 100, 500, 800, 1000, 1500, 2000, 3000];
for kk = 1:length(krho_tests)
    kr = krho_tests(kk) + 0.001j;
    kz_t = sqrt(k^2 - kr^2);
    if imag(kz_t)>0, kz_t = complex(real(kz_t),-abs(imag(kz_t))); end
    if real(kz_t)<0, kz_t = complex(abs(real(kz_t)),imag(kz_t)); end
    
    smgf.SetRadialWaveNumber(kr);
    smgf.ComputeSpectralMGF();
    K_dir = smgf.K(comp) * kz_t;
    
    K_dcim = 0;
    for jj = 1:length(all_a)
        K_dcim = K_dcim + all_a(jj) * exp(-all_alpha(jj) * kz_t);
    end
    
    err = abs(K_dcim - K_dir)/abs(K_dir)*100;
    fprintf('krho=%5.0f: |K_dir|=%.3e, |K_dcim|=%.3e, err=%.1f%%\n', ...
        real(kr), abs(K_dir), abs(K_dcim), err);
end

%% 最终对比
G_dcim = 0;
for ii = 1:length(all_a)
    Rc = sqrt(rho_test^2 - all_alpha(ii)^2);
    G_dcim = G_dcim + all_a(ii) * exp(-1j*k*Rc) / Rc;
end
G_dcim = G_dcim * 1j / (2*pi);

% INTEGRATE 基准
smgf2 = SpectralMGF();
smgf2.Initialize(lm, freq, [1,1,1,1,1], true, false);
smgf2.SetLayers(i_src, m_obs);
smgf2.SetSourcePoint(z_src);
smgf2.SetObservationPoint(z_obs);
a_sw = 1.2 * abs(lm.k_max);

G_int = SommerfeldIntegrator2.IntegrateSpectralNearField(smgf2, rho_test, comp, 0, a_sw, false) ...
      + SommerfeldIntegrator2.IntegrateSpectralFarField(smgf2, rho_test, comp, 0, a_sw, false);

% 准静态
qmgf = QuasistaticMGF2();
qmgf.Initialize(lm, freq, [1,1,1,1,1], false);
qmgf.SetLayers(i_src, m_obs);
K_quasi = qmgf.ComputeQMGF_Spatial(z_obs, z_src, rho_test);
G_int_total = G_int + K_quasi(comp);

fprintf('\n=== 结果 ===\n');
fprintf('G_dcim      = %.6e %+.6ei\n', real(G_dcim), imag(G_dcim));
fprintf('G_int_total = %.6e %+.6ei\n', real(G_int_total), imag(G_int_total));
fprintf('误差 = %.2f%%\n', abs(G_dcim - G_int_total)/abs(G_int_total)*100);
