%% 对比不同 GPOF 实现策略
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

t2 = linspace(1e-2, T02, N2);
kz2 = k*(-1j*t2+(1-t2/T02));
krho2 = sqrt(k^2-kz2.^2); krho2 = abs(real(krho2))+1j*abs(imag(krho2));

% 采样路径2（直接做单层 DCIM 简化）
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
fprintf('G_ref = %.4e %+.4ei (abs=%.4f)\n\n', real(G_ref), imag(G_ref), abs(G_ref));

%% 方法A: 当前 GPOF (full SVD + full Z)
fprintf('=== 方法A: 当前 GPOF ===\n');
[at_A, a_A, res_A] = GPOF(K2(:), t2(:), 50, 1e-4, 1e-16);
fprintf('%d exp, res=%.4e\n', length(at_A), res_A);

al_kz_A = at_A.*T02./((1+1j*T02)*k);
a_kz_A = a_A.*exp(k*al_kz_A);
G_A = compute_G(a_kz_A, al_kz_A, k, rho_test);
fprintf('G = %.4e %+.4ei (abs=%.4f, err=%.2f%%)\n', real(G_A), imag(G_A), abs(G_A), abs(G_A-G_ref)/abs(G_ref)*100);

%% 方法B: pinv 代替 SVD 截断
fprintf('\n=== 方法B: pinv-based GPOF ===\n');
y = K2(:); t = t2(:); N = length(y); L = 50;
dt = t(2)-t(1);
m_rows = N-L; n_cols = L;

Y1 = zeros(m_rows, n_cols);
Y2_mat = zeros(m_rows, n_cols);
for ii=1:n_cols
    for jj=1:m_rows
        Y1(jj,ii) = y(jj+ii-1);
        Y2_mat(jj,ii) = y(jj+ii);
    end
end

% 用 pinv 方法
Z_pinv = pinv(Y1) * Y2_mat;
w_B = eig(Z_pinv);
[~,idx] = sort(abs(w_B),'descend');
w_B = w_B(idx);
% 过滤极小特征值
mask = abs(w_B) > 1e-4 * abs(w_B(1));
w_B = w_B(mask);

Y3 = zeros(N, length(w_B));
for ii=1:length(w_B), Y3(:,ii) = w_B(ii).^(0:N-1)'; end
a_B = Y3 \ y;
at_B = log(w_B)/dt;

res_B = norm(y - Y3*a_B)/norm(y);
fprintf('%d exp, res=%.4e\n', length(at_B), res_B);

al_kz_B = at_B.*T02./((1+1j*T02)*k);
a_kz_B = a_B.*exp(k*al_kz_B);
G_B = compute_G(a_kz_B, al_kz_B, k, rho_test);
fprintf('G = %.4e %+.4ei (abs=%.4f, err=%.2f%%)\n', real(G_B), imag(G_B), abs(G_B), abs(G_B-G_ref)/abs(G_ref)*100);

%% 方法C: TLS-ESPRIT 方式 (不同的矩阵构造)
fprintf('\n=== 方法C: 广义特征值 eig(Y2,Y1) ===\n');
% 直接用广义特征值
w_C = eig(Y2_mat, Y1);
w_C = w_C(isfinite(w_C));
[~,idx] = sort(abs(w_C),'descend');
w_C = w_C(idx);

% 用 tol_svd 对应的方式截断：只保留 |w| > threshold
% 使用 SVD 确定有效秩
sv = svd(Y1);
nst = sum(sv > 1e-4*sv(1));
w_C = w_C(1:min(nst, length(w_C)));

Y3c = zeros(N, length(w_C));
for ii=1:length(w_C), Y3c(:,ii) = w_C(ii).^(0:N-1)'; end
a_C = Y3c \ y;
at_C = log(w_C)/dt;

res_C = norm(y - Y3c*a_C)/norm(y);
fprintf('%d exp, res=%.4e\n', length(at_C), res_C);

al_kz_C = at_C.*T02./((1+1j*T02)*k);
a_kz_C = a_C.*exp(k*al_kz_C);
G_C = compute_G(a_kz_C, al_kz_C, k, rho_test);
fprintf('G = %.4e %+.4ei (abs=%.4f, err=%.2f%%)\n', real(G_C), imag(G_C), abs(G_C), abs(G_C-G_ref)/abs(G_ref)*100);

%% 方法D: 遍历不同指数数量，找最优
fprintf('\n=== 方法D: 扫描指数数量（单层DCIM简化） ===\n');
[U,S,V] = svd(Y1);
sv = diag(S);
ns = min(m_rows, n_cols);

for M_test = 2:min(15, ns)
    U_M = U(:,1:M_test); V_M = V(:,1:M_test); S_M = S(1:M_test,1:M_test);
    Z_M = S_M \ (U_M' * Y2_mat * V_M);
    w_M = eig(Z_M);
    [~,idx] = sort(abs(w_M),'descend');
    w_M = w_M(idx);
    
    Y3m = zeros(N, M_test);
    for ii=1:M_test, Y3m(:,ii) = w_M(ii).^(0:N-1)'; end
    a_M = Y3m \ y;
    at_M = log(w_M)/dt;
    res_M = norm(y - Y3m*a_M)/norm(y);
    
    al_kz_M = at_M.*T02./((1+1j*T02)*k);
    a_kz_M = a_M.*exp(k*al_kz_M);
    G_M = compute_G(a_kz_M, al_kz_M, k, rho_test);
    err_M = abs(G_M-G_ref)/abs(G_ref)*100;
    fprintf('M=%2d: res=%.3e, max|a_kz|=%.2e, abs(G)=%.2f, err=%6.2f%%\n', ...
        M_test, res_M, max(abs(a_kz_M)), abs(G_M), err_M);
end

function G = compute_G(aa, aal, k, rho)
    G = 0;
    for ii=1:length(aa)
        Rc = sqrt(rho^2-aal(ii)^2);
        G = G + aa(ii)*exp(-1j*k*Rc)/Rc;
    end
    G = G*1j/(2*pi);
end
