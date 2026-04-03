%% 诊断DCIM问题
clear; clc;
addpath(pwd);

%% 设置
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

fprintf('k_min = %.6f + %.6fi\n', real(lm.k_min), imag(lm.k_min));
fprintf('k_max = %.6f + %.6fi\n', real(lm.k_max), imag(lm.k_max));
fprintf('k_top = %.6f + %.6fi\n', real(lm.k_top), imag(lm.k_top));
fprintf('k_bot = %.6f + %.6fi\n', real(lm.k_bot), imag(lm.k_bot));
for ii = 1:length(lm.layers)
    fprintf('层%d: k = %.6f, epsr = %.1f\n', ii, lm.layers(ii).k, lm.layers(ii).epsr);
end

%% DCIM 参数
k = lm.k_min;
epsr_list = [lm.layers.epsr];
epsr_max = max(epsr_list);
T02 = 5.0 * sqrt(epsr_max);
T01 = 200.0;
N1 = 100;
N2 = 100;

fprintf('\n=== DCIM 参数 ===\n');
fprintf('k = %.6f + %.6fi\n', real(k), imag(k));
fprintf('epsr_max = %.2f\n', epsr_max);
fprintf('T02 = %.4f\n', T02);
fprintf('T01 = %.4f\n', T01);

%% 生成路径1
t1 = linspace(0, T01, N1);
path1_kz = -1j * k * (T02 + t1);
path1_krho = sqrt(k^2 - path1_kz.^2);
path1_krho = abs(real(path1_krho)) + 1j * abs(imag(path1_krho));

fprintf('\n=== 路径1 ===\n');
fprintf('path1_kz(1)  = %.6e + %.6ei\n', real(path1_kz(1)), imag(path1_kz(1)));
fprintf('path1_kz(end)= %.6e + %.6ei\n', real(path1_kz(end)), imag(path1_kz(end)));
fprintf('path1_krho(1)  = %.6e + %.6ei\n', real(path1_krho(1)), imag(path1_krho(1)));
fprintf('path1_krho(end)= %.6e + %.6ei\n', real(path1_krho(end)), imag(path1_krho(end)));

%% 生成路径2
t2 = linspace(1e-2, T02, N2);
path2_kz = k * (-1j * t2 + (1 - t2/T02));
path2_krho = sqrt(k^2 - path2_kz.^2);
path2_krho = abs(real(path2_krho)) + 1j * abs(imag(path2_krho));

fprintf('\n=== 路径2 ===\n');
fprintf('path2_kz(1)  = %.6e + %.6ei\n', real(path2_kz(1)), imag(path2_kz(1)));
fprintf('path2_kz(end)= %.6e + %.6ei\n', real(path2_kz(end)), imag(path2_kz(end)));
fprintf('path2_krho(1)  = %.6e + %.6ei\n', real(path2_krho(1)), imag(path2_krho(1)));
fprintf('path2_krho(end)= %.6e + %.6ei\n', real(path2_krho(end)), imag(path2_krho(end)));

%% 采样频谱域格林函数（路径1，分量1=Gxx）
smgf = SpectralMGF();
smgf.Initialize(lm, freq, [1,1,1,1,1], true, false);
smgf.SetLayers(i_src, m_obs);
smgf.SetSourcePoint(z_src);
smgf.SetObservationPoint(z_obs);

comp_idx = 1;
K_path1 = zeros(1, N1);
for ii = 1:N1
    smgf.SetRadialWaveNumber(path1_krho(ii));
    smgf.ComputeSpectralMGF();
    K_path1(ii) = smgf.K(comp_idx) * path1_kz(ii);
end

fprintf('\n=== 路径1采样结果 (分量 %d) ===\n', comp_idx);
fprintf('K_path1(1) = %.6e + %.6ei\n', real(K_path1(1)), imag(K_path1(1)));
fprintf('K_path1(50) = %.6e + %.6ei\n', real(K_path1(50)), imag(K_path1(50)));
fprintf('K_path1(end) = %.6e + %.6ei\n', real(K_path1(end)), imag(K_path1(end)));
fprintf('max(abs(K_path1)) = %.6e\n', max(abs(K_path1)));
fprintf('min(abs(K_path1)) = %.6e\n', min(abs(K_path1)));
fprintf('norm(K_path1) = %.6e\n', norm(K_path1));

%% GPOF 拟合路径1
dt1 = t1(2) - t1(1);
L1 = 50;
fprintf('\ndt1 = %.6e\n', dt1);

[alpha_t1, a1, res1] = GPOF(K_path1(:), t1(:), L1, 1e-4, 1e-16);
fprintf('\n=== GPOF Level 1 结果 ===\n');
fprintf('GPOF 找到 %d 个指数\n', length(alpha_t1));
fprintf('残差 = %.6e\n', res1);
for jj = 1:min(5, length(alpha_t1))
    fprintf('  alpha_t(%d) = %.6e + %.6ei, a(%d) = %.6e + %.6ei\n', ...
        jj, real(alpha_t1(jj)), imag(alpha_t1(jj)), jj, real(a1(jj)), imag(a1(jj)));
end

%% 转换到kz空间
a1_kz = a1 .* exp(-T02 * alpha_t1);
alpha1_kz = alpha_t1 ./ (1j * k);

fprintf('\n=== Level 1 转换后 ===\n');
for jj = 1:min(5, length(a1_kz))
    fprintf('  a_kz(%d) = %.6e + %.6ei, alpha_kz(%d) = %.6e + %.6ei\n', ...
        jj, real(a1_kz(jj)), imag(a1_kz(jj)), jj, real(alpha1_kz(jj)), imag(alpha1_kz(jj)));
end

%% 验证Level 1拟合质量
K_path1_fit = zeros(1, N1);
for ii = 1:N1
    for jj = 1:length(a1)
        K_path1_fit(ii) = K_path1_fit(ii) + a1(jj) * exp(alpha_t1(jj) * t1(ii));
    end
end
fit_err1 = norm(K_path1 - K_path1_fit) / norm(K_path1);
fprintf('Level 1 拟合误差 = %.6e\n', fit_err1);

%% 采样路径2并减去Level 1
K_path2 = zeros(1, N2);
for ii = 1:N2
    smgf.SetRadialWaveNumber(path2_krho(ii));
    smgf.ComputeSpectralMGF();
    K_path2(ii) = smgf.K(comp_idx) * path2_kz(ii);
end

fprintf('\n=== 路径2采样结果 (分量 %d) ===\n', comp_idx);
fprintf('K_path2(1) = %.6e + %.6ei\n', real(K_path2(1)), imag(K_path2(1)));
fprintf('K_path2(50) = %.6e + %.6ei\n', real(K_path2(50)), imag(K_path2(50)));
fprintf('K_path2(end) = %.6e + %.6ei\n', real(K_path2(end)), imag(K_path2(end)));
fprintf('max(abs(K_path2)) = %.6e\n', max(abs(K_path2)));
fprintf('norm(K_path2) = %.6e\n', norm(K_path2));

% 减去Level 1
K_path2_mod = K_path2;
for ii = 1:N2
    kz = path2_kz(ii);
    for jj = 1:length(a1_kz)
        K_path2_mod(ii) = K_path2_mod(ii) - a1_kz(jj) * exp(-alpha1_kz(jj) * kz);
    end
end

fprintf('\n减去Level 1后:\n');
fprintf('max(abs(K_path2_mod)) = %.6e\n', max(abs(K_path2_mod)));
fprintf('norm(K_path2_mod) = %.6e\n', norm(K_path2_mod));
fprintf('norm(K_path2_mod)/norm(K_path2) = %.6e\n', norm(K_path2_mod)/norm(K_path2));

%% GPOF 拟合路径2
dt2 = t2(2) - t2(1);
[alpha_t2, a2, res2] = GPOF(K_path2_mod(:), t2(:), 50, 1e-4, 1e-16);
fprintf('\n=== GPOF Level 2 结果 ===\n');
fprintf('GPOF 找到 %d 个指数\n', length(alpha_t2));
fprintf('残差 = %.6e\n', res2);

% 转换
alpha2_kz = alpha_t2 .* T02 ./ ((1 + 1j*T02) * k);
a2_kz = a2 .* exp(k * alpha2_kz);

fprintf('\n=== Level 2 转换后 ===\n');
for jj = 1:min(5, length(a2_kz))
    fprintf('  a_kz(%d) = %.6e + %.6ei, alpha_kz(%d) = %.6e + %.6ei\n', ...
        jj, real(a2_kz(jj)), imag(a2_kz(jj)), jj, real(alpha2_kz(jj)), imag(alpha2_kz(jj)));
end

%% 合并并计算 G
all_a = [a1_kz(:); a2_kz(:)];
all_alpha = [alpha1_kz(:); alpha2_kz(:)];
rho_test = 1e-6;

% Sommerfeld Identity (S0 for comp 1 = Gxx)
G_dcim_core = 0;
for ii = 1:length(all_a)
    Rc = sqrt(rho_test^2 - all_alpha(ii)^2);
    G_dcim_core = G_dcim_core + all_a(ii) * exp(-1j*k*Rc) / Rc;
end
G_dcim_core = G_dcim_core * 1j / (2*pi);

% 加回准静态项
qmgf = QuasistaticMGF2();
qmgf.Initialize(lm, freq, [1,1,1,1,1], false);
qmgf.SetLayers(i_src, m_obs);
K_quasi = qmgf.ComputeQMGF_Spatial(z_obs, z_src, rho_test);

G_dcim_total = G_dcim_core + K_quasi(comp_idx);

fprintf('\n=== 最终结果 (分量 %d) ===\n', comp_idx);
fprintf('G_dcim_core = %.6e + %.6ei\n', real(G_dcim_core), imag(G_dcim_core));
fprintf('K_quasi(%d) = %.6e + %.6ei\n', comp_idx, real(K_quasi(comp_idx)), imag(K_quasi(comp_idx)));
fprintf('G_dcim_total = %.6e + %.6ei\n', real(G_dcim_total), imag(G_dcim_total));

%% 基准
smgf2 = SpectralMGF();
smgf2.Initialize(lm, freq, [1,1,1,1,1], true, false);
smgf2.SetLayers(i_src, m_obs);
smgf2.SetSourcePoint(z_src);
smgf2.SetObservationPoint(z_obs);

a_switch = 1.2 * abs(lm.k_max);
G_int = SommerfeldIntegrator2.IntegrateSpectralNearField(smgf2, rho_test, comp_idx, 0, a_switch, false) ...
      + SommerfeldIntegrator2.IntegrateSpectralFarField(smgf2, rho_test, comp_idx, 0, a_switch, false);
G_int_total = G_int + K_quasi(comp_idx);

fprintf('\nG_integrate_core = %.6e + %.6ei\n', real(G_int), imag(G_int));
fprintf('G_integrate_total = %.6e + %.6ei\n', real(G_int_total), imag(G_int_total));
fprintf('误差 = %.2f%%\n', abs(G_dcim_total - G_int_total)/abs(G_int_total)*100);
