%% 简化诊断 - 对比 DCIM 和 INTEGRATE 在相同 rho 下的表现
clear; clc;
addpath(pwd);

freq = 30e9;
unit = 1e-3;
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

fprintf('k_min = %.6f + %.6fi (abs=%.6f)\n', real(lm.k_min), imag(lm.k_min), abs(lm.k_min));
fprintf('k_max = %.6f + %.6fi (abs=%.6f)\n', real(lm.k_max), imag(lm.k_max), abs(lm.k_max));

%% DCIM 手动实现（简化，方便诊断）
k = lm.k_min;
epsr_list = [lm.layers.epsr];
epsr_max = max(epsr_list);
T02 = 5.0 * sqrt(epsr_max);
T01 = 200.0;
N1 = 100; N2 = 100;

fprintf('\nk = %.6f, T02 = %.4f, T01 = %.4f\n', k, T02, T01);

% 创建 SpectralMGF
smgf = SpectralMGF();
smgf.Initialize(lm, freq, [1,1,1,1,1], true, false);
smgf.SetLayers(i_src, m_obs);
smgf.SetSourcePoint(z_src);
smgf.SetObservationPoint(z_obs);

% 路径1
t1 = linspace(0, T01, N1);
kz1 = -1j * k * (T02 + t1);
krho1 = sqrt(k^2 - kz1.^2);
krho1 = abs(real(krho1)) + 1j * abs(imag(krho1));

% 路径2
t2 = linspace(1e-2, T02, N2);
kz2 = k * (-1j * t2 + (1 - t2/T02));
krho2 = sqrt(k^2 - kz2.^2);
krho2 = abs(real(krho2)) + 1j * abs(imag(krho2));

% 采样所有分量 (只看分量1=Gxx)
comp = 1;

% 路径1采样
K1 = zeros(1, N1);
for ii = 1:N1
    smgf.SetRadialWaveNumber(krho1(ii));
    smgf.ComputeSpectralMGF();
    K1(ii) = smgf.K(comp) * kz1(ii);
end
fprintf('\n路径1: max|K|=%.4e, min|K|=%.4e\n', max(abs(K1)), min(abs(K1)));

% GPOF Level 1
[alpha_t1, a1, res1] = GPOF(K1(:), t1(:), 50, 1e-4, 1e-16);
fprintf('GPOF L1: %d 指数, 残差=%.4e\n', length(alpha_t1), res1);

% t->kz 转换 (Level 1)
a1_kz = a1 .* exp(-T02 * alpha_t1);
alpha1_kz = alpha_t1 ./ (1j * k);
fprintf('L1 转换后: max|a|=%.4e\n', max(abs(a1_kz)));

% 路径2采样
K2 = zeros(1, N2);
for ii = 1:N2
    smgf.SetRadialWaveNumber(krho2(ii));
    smgf.ComputeSpectralMGF();
    K2(ii) = smgf.K(comp) * kz2(ii);
end
fprintf('路径2: max|K|=%.4e, min|K|=%.4e\n', max(abs(K2)), min(abs(K2)));

% 减去 Level 1 贡献
K2_mod = K2;
for ii = 1:N2
    for jj = 1:length(a1_kz)
        K2_mod(ii) = K2_mod(ii) - a1_kz(jj) * exp(-alpha1_kz(jj) * kz2(ii));
    end
end
fprintf('减去L1后: max|K2_mod|=%.4e (%.1f%% of original)\n', ...
    max(abs(K2_mod)), max(abs(K2_mod))/max(abs(K2))*100);

% GPOF Level 2
[alpha_t2, a2, res2] = GPOF(K2_mod(:), t2(:), 50, 1e-4, 1e-16);
fprintf('GPOF L2: %d 指数, 残差=%.4e\n', length(alpha_t2), res2);

% t->kz 转换 (Level 2)
alpha2_kz = alpha_t2 .* T02 ./ ((1 + 1j*T02) * k);
a2_kz = a2 .* exp(k * alpha2_kz);
fprintf('L2 转换后: max|a|=%.4e\n', max(abs(a2_kz)));

% 合并
all_a = [a1_kz(:); a2_kz(:)];
all_alpha = [alpha1_kz(:); alpha2_kz(:)];

%% 验证: 在路径2上重建并与原始数据对比
K2_recon = zeros(1, N2);
for ii = 1:N2
    kz = kz2(ii);
    for jj = 1:length(all_a)
        K2_recon(ii) = K2_recon(ii) + all_a(jj) * exp(-all_alpha(jj) * kz);
    end
end
fit_err = norm(K2 - K2_recon) / norm(K2);
fprintf('\n路径2重建误差 = %.4e\n', fit_err);

%% 计算空间域格林函数
rho_test = 1e-6;

% DCIM Sommerfeld identity
G_dcim = 0;
for ii = 1:length(all_a)
    Rc = sqrt(rho_test^2 - all_alpha(ii)^2);
    G_dcim = G_dcim + all_a(ii) * exp(-1j*k*Rc) / Rc;
end
G_dcim = G_dcim * 1j / (2*pi);

% 准静态项
qmgf = QuasistaticMGF2();
qmgf.Initialize(lm, freq, [1,1,1,1,1], false);
qmgf.SetLayers(i_src, m_obs);
K_quasi = qmgf.ComputeQMGF_Spatial(z_obs, z_src, rho_test);

G_dcim_total = G_dcim + K_quasi(comp);

% INTEGRATE
smgf2 = SpectralMGF();
smgf2.Initialize(lm, freq, [1,1,1,1,1], true, false);
smgf2.SetLayers(i_src, m_obs);
smgf2.SetSourcePoint(z_src);
smgf2.SetObservationPoint(z_obs);
a_switch = 1.2 * abs(lm.k_max);

G_int = SommerfeldIntegrator2.IntegrateSpectralNearField(smgf2, rho_test, comp, 0, a_switch, false) ...
      + SommerfeldIntegrator2.IntegrateSpectralFarField(smgf2, rho_test, comp, 0, a_switch, false);
G_int_total = G_int + K_quasi(comp);

fprintf('\n=== 结果对比 ===\n');
fprintf('G_dcim_core  = %.6e %+.6ei\n', real(G_dcim), imag(G_dcim));
fprintf('G_int_core   = %.6e %+.6ei\n', real(G_int), imag(G_int));
fprintf('K_quasi      = %.6e %+.6ei\n', real(K_quasi(comp)), imag(K_quasi(comp)));
fprintf('G_dcim_total = %.6e %+.6ei\n', real(G_dcim_total), imag(G_dcim_total));
fprintf('G_int_total  = %.6e %+.6ei\n', real(G_int_total), imag(G_int_total));
fprintf('误差 = %.2f%%\n', abs(G_dcim_total - G_int_total)/abs(G_int_total)*100);

%% 额外诊断: 在实轴上验证 DCIM 拟合
fprintf('\n=== 实轴验证 ===\n');
krho_test_vals = [1, 100, 500, 1000, 2000] + 0.001j;
for kk = 1:length(krho_test_vals)
    kr = krho_test_vals(kk);
    kz_test = sqrt(k^2 - kr^2);
    if imag(kz_test) > 0
        kz_test = complex(real(kz_test), -abs(imag(kz_test)));
    end
    if real(kz_test) < 0
        kz_test = complex(abs(real(kz_test)), imag(kz_test));
    end
    
    % 直接计算
    smgf.SetRadialWaveNumber(kr);
    smgf.ComputeSpectralMGF();
    K_direct = smgf.K(comp) * kz_test;
    
    % DCIM 重建
    K_dcim_val = 0;
    for jj = 1:length(all_a)
        K_dcim_val = K_dcim_val + all_a(jj) * exp(-all_alpha(jj) * kz_test);
    end
    
    fprintf('krho=%.0f: K_direct=%.4e%+.4ei, K_dcim=%.4e%+.4ei, err=%.1f%%\n', ...
        real(kr), real(K_direct), imag(K_direct), real(K_dcim_val), imag(K_dcim_val), ...
        abs(K_dcim_val - K_direct)/abs(K_direct)*100);
end
