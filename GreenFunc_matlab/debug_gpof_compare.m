%% GPOF 数值对比：MATLAB vs C++ strata
% 输出格式与 C++ debug 输出完全对应
clear classes; clear; clc;

%% 设置（与 main_test_dcim_standard.m 和 testDCIM.cpp 一致）
freq = 30e9;
eps0 = Constants.eps0; mu0 = Constants.mu0;

lm = LayerManager();
lm.SetTopHalfspace(1.0, 1.0, false);
lm.SetBottomHalfspace(1.0, 1.0, true);
lm.AddLayer(0, 0.4e-3, 2.1, 1, 0, 0);
lm.AddLayer(0.4e-3, 0.8e-3, 12.5, 1, 0, 0);
lm.AddLayer(0.8e-3, 1.4e-3, 9.8, 1, 0, 0);
lm.AddLayer(1.4e-3, 1.9e-3, 8.6, 1, 0, 0);
lm.ProcessLayers(freq);

z_src = 0.4e-3; z_obs = 1.4e-3;
i_src = 2; m_obs = 0;
omega = 2*pi*freq;

%% DCIM 参数（与 strata 默认一致）
k = lm.k_min;
epsr_list = [lm.layers.epsr];
epsr_max = max(epsr_list);
T02 = 5.0 * sqrt(epsr_max);
T01 = 200.0;
N1 = 100; N2 = 100;
tol_svd = 1e-4; tol_eig = 1e-16;

%% 生成路径
t1 = linspace(0, T01, N1);
kz1 = -1j * k * (T02 + t1);
krho1 = sqrt(k^2 - kz1.^2);
krho1 = abs(real(krho1)) + 1j*abs(imag(krho1));

t2 = linspace(1e-2, T02, N2);
kz2 = k * (-1j * t2 + (1 - t2/T02));
krho2 = sqrt(k^2 - kz2.^2);
krho2 = abs(real(krho2)) + 1j*abs(imag(krho2));

dt1 = t1(2) - t1(1);
dt2 = t2(2) - t2(1);

%% 输出路径参数（对应 C++ "DCIM DEBUG"）
fprintf('==== DCIM DEBUG (comp 0) ====\n');
fprintf('k = (%.15e, %.15e)\n', real(k), imag(k));
fprintf('T02 = %.15e, T01 = %.15e\n', T02, T01);
fprintf('N1 = %d, N2 = %d\n', N1, N2);
fprintf('dt1 = %.15e, dt2 = %.15e\n', dt1, dt2);

%% 采样谱域 MGF（comp 0 = Gxx, extract_quasistatic=false）
smgf = SpectralMGF();
smgf.Initialize(lm, freq, [1,1,1,1,1], false, false);
smgf.SetLayers(i_src, m_obs);
smgf.SetSourcePoint(z_src);
smgf.SetObservationPoint(z_obs);

K_path1 = zeros(1, N1);
for mm = 1:N1
    smgf.SetRadialWaveNumber(krho1(mm));
    smgf.ComputeSpectralMGF();
    K_path1(mm) = smgf.K(1) * kz1(mm);
end

K_path2 = zeros(1, N2);
for mm = 1:N2
    smgf.SetRadialWaveNumber(krho2(mm));
    smgf.ComputeSpectralMGF();
    K_path2(mm) = smgf.K(1) * kz2(mm);
end

%% 输出路径2采样值（对应 C++ "K_path2 before L1 subtract"）
fprintf('\n--- K_path2 (before L1 subtract) first 5 ---\n');
for pp = 1:5
    fprintf('  K_path2[%d] = (%.15e, %.15e)\n', pp-1, real(K_path2(pp)), imag(K_path2(pp)));
end
fprintf('  K_path2[%d] = (%.15e, %.15e)\n', N2-1, real(K_path2(end)), imag(K_path2(end)));

%% ============ GPOF Level 1 ============
fprintf('\n--- GPOF call #1 ---\n');
y1 = K_path1(:);
[alpha_t1, a_t1, svd_info1] = GPOF_debug(y1, t1(:), N1/2, tol_svd, tol_eig);

% t->kz 转换
if ~isempty(a_t1)
    a1_kz = a_t1 .* exp(-T02 * alpha_t1);
    alpha1_kz = alpha_t1 ./ (1j * k);
else
    a1_kz = []; alpha1_kz = [];
end

fprintf('\n--- Level 1 after t->kz conversion ---\n');
for jj = 1:length(a1_kz)
    fprintf('  a_kz[%d]=(%.15e,%.15e)  alpha_kz[%d]=(%.15e,%.15e)\n', ...
        jj-1, real(a1_kz(jj)), imag(a1_kz(jj)), jj-1, real(alpha1_kz(jj)), imag(alpha1_kz(jj)));
end

%% 从路径2减去 Level 1
K2_mod = K_path2;
for mm = 1:N2
    for jj = 1:length(a1_kz)
        K2_mod(mm) = K2_mod(mm) - a1_kz(jj) * exp(-alpha1_kz(jj) * kz2(mm));
    end
end

%% ============ GPOF Level 2 ============
fprintf('\n--- GPOF call #2 ---\n');
y2 = K2_mod(:);
[alpha_t2, a_t2, svd_info2] = GPOF_debug(y2, t2(:), N2/2, tol_svd, tol_eig);

% t->kz 转换
if ~isempty(a_t2)
    alpha2_kz = alpha_t2 .* T02 ./ ((1 + 1j*T02) * k);
    a2_kz = a_t2 .* exp(k * alpha2_kz);
else
    a2_kz = []; alpha2_kz = [];
end

fprintf('\n--- Level 2 after t->kz conversion ---\n');
for jj = 1:length(a2_kz)
    fprintf('  a_kz[%d]=(%.15e,%.15e)  alpha_kz[%d]=(%.15e,%.15e)\n', ...
        jj-1, real(a2_kz(jj)), imag(a2_kz(jj)), jj-1, real(alpha2_kz(jj)), imag(alpha2_kz(jj)));
end

%% ======================================
%% GPOF_debug: 带完整中间值输出的 GPOF
%% ======================================
function [exponents, coefficients, info] = GPOF_debug(y, t, L, tol_svd, tol_eig)

    N = length(y);
    dt = t(2) - t(1);
    m = N - L;  % rows
    n = L;      % cols
    
    % 构造 Hankel 矩阵 Y1, Y2
    Y1 = zeros(m, n);
    Y2 = zeros(m, n);
    for ii = 1:n
        for jj = 1:m
            Y1(jj, ii) = y(jj + ii - 1);
            Y2(jj, ii) = y(jj + ii);
        end
    end
    
    % Full SVD (与 C++ LAPACK 'A','A' 一致)
    [U, S_mat, V] = svd(Y1);
    s = diag(S_mat);
    ns = min(m, n);
    s = s(1:ns);
    U_ns = U(:, 1:ns);
    V_ns = V(:, 1:ns);
    
    fprintf('N=%d, L=%d, m=%d, n=%d, ns=%d\n', N, L, m, n, ns);
    fprintf('Singular values (first 10):\n');
    for pp = 1:min(10, ns)
        fprintf('  s[%d] = %.15e\n', pp-1, s(pp));
    end
    
    % SVD 截断 (与 C++ 完全对应)
    nst = ns;  % 默认保留全部
    for pp = 2:ns
        if s(pp) < tol_svd * s(1)
            nst = pp - 1;
            break;
        end
    end
    fprintf('nst (SVD cutoff) = %d\n', nst);
    
    % D_inv: ns x ns 对角矩阵，nst 之后设为 0
    D_inv = zeros(ns, ns);
    for pp = 1:nst
        D_inv(pp, pp) = 1.0 / s(pp);
    end
    
    % Z = D_inv * U_ns' * Y2 * V_ns (ns x ns)
    Z = D_inv * (U_ns' * Y2 * V_ns);
    
    % 特征值
    w = eig(Z);
    
    % 按绝对值降序排列（C++ zgeev 不保证排序，但我们需要对齐对比）
    [~, idx] = sort(abs(w), 'descend');
    w = w(idx);
    
    fprintf('Eigenvalues of Z (first 10):\n');
    for pp = 1:min(10, ns)
        fprintf('  w[%d] = (%.15e, %.15e), |w|=%.15e\n', pp-1, real(w(pp)), imag(w(pp)), abs(w(pp)));
    end
    
    % 特征值截断
    nwt = ns;
    for pp = 2:ns
        if abs(w(pp)) < tol_eig * abs(w(1))
            nwt = pp - 1;
            break;
        end
    end
    fprintf('nwt (eig cutoff) = %d\n', nwt);
    
    w = w(1:nwt);
    
    % 构造 Vandermonde 矩阵并求解系数
    Y3 = zeros(N, nwt);
    for ii = 1:nwt
        for jj = 1:N
            Y3(jj, ii) = w(ii)^(jj-1);
        end
    end
    
    b = Y3 \ y;
    
    fprintf('Final: %d exponents\n', nwt);
    for pp = 1:nwt
        fprintf('  a[%d]=(%.15e,%.15e)  alpha[%d]=(%.15e,%.15e)\n', ...
            pp-1, real(b(pp)), imag(b(pp)), pp-1, real(log(w(pp))/dt), imag(log(w(pp))/dt));
    end
    
    exponents = log(w) / dt;
    coefficients = b;
    info.s = s;
    info.nst = nst;
    info.nwt = nwt;
    info.w = w;
end
