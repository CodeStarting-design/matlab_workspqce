%% 测试DCIM方法（标准配置）
clear; clc;
addpath(pwd);

fprintf('========================================\n');
fprintf('DCIM测试（标准配置）\n');
fprintf('========================================\n\n');

%% 标准配置
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
rho_test = 1e-6;

i_src = lm.FindLayer(z_src);
m_obs = lm.FindLayer(z_obs);

fprintf('配置: f=%.0f GHz, z_src=%.2f mm (idx=%d), z_obs=%.2f mm (idx=%d)\n\n', ...
    freq/1e9, z_src/unit, i_src, z_obs/unit, m_obs);

%% 方法1: INTEGRATE（基准）
fprintf('方法1: INTEGRATE（基准）\n');
s1 = struct();
s1.method = 'INTEGRATE';
s1.components = [1,1,1,1,1];
s1.extract_quasistatic = true;

mgf1 = MGF3_Integrated();
mgf1.Initialize(freq, lm, s1);
mgf1.SetLayers(i_src, m_obs);

tic;
[G1, ~] = mgf1.ComputeMGF(rho_test, 0, z_obs, z_src);
t1 = toc;
G_integrate = G1(1);
fprintf('G_integrate = %+.4e %+.4ej (%.3f秒)\n\n', real(G_integrate), imag(G_integrate), t1);

%% 方法2: DCIM
fprintf('方法2: DCIM\n');
s2 = struct();
s2.method = 'DCIM';
s2.components = [1,1,1,1,1];
s2.extract_quasistatic = true;

mgf2 = MGF3_Integrated();
mgf2.Initialize(freq, lm, s2);

tic;
mgf2.SetLayers(i_src, m_obs, z_src, z_obs);
t2_setup = toc;

tic;
[G2, ~] = mgf2.ComputeMGF(rho_test, 0, z_obs, z_src);
t2_compute = toc;

G_dcim = G2(1);
fprintf('G_dcim = %+.4e %+.4ej (%.3f秒设置 + %.3f秒计算)\n\n', ...
    real(G_dcim), imag(G_dcim), t2_setup, t2_compute);

%% 对比
err = abs(G_dcim - G_integrate) / abs(G_integrate) * 100;
fprintf('误差: %.2f%%\n', err);
