%% 两层 DCIM + 自动选择最优 M（通过 max|a_kz| 限制）
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

%% 两层 DCIM，Level 2 用自适应 M 选择
% Level 1 GPOF
[at1, a1, ~] = GPOF(K1(:), t1(:), 50, 1e-4, 1e-16);
a1k = a1.*exp(-T02*at1); al1k = at1./(1j*k);
fprintf('L1: %d exp\n', length(at1));

% 减去 Level 1
K2m = K2;
for ii=1:N2
    for jj=1:length(a1k)
        K2m(ii) = K2m(ii)-a1k(jj)*exp(-al1k(jj)*kz2(ii));
    end
end

% Level 2: 扫描 M 选择最优
y2 = K2m(:); t2v = t2(:);
N = length(y2); L = 50;
dt2 = t2v(2)-t2v(1);
m_rows = N-L; n_cols = L;

Y1mat = zeros(m_rows, n_cols);
Y2mat = zeros(m_rows, n_cols);
for ii=1:n_cols
    for jj=1:m_rows
        Y1mat(jj,ii) = y2(jj+ii-1);
        Y2mat(jj,ii) = y2(jj+ii);
    end
end

[U,S,V] = svd(Y1mat);
sv = diag(S);
ns = min(m_rows, n_cols);

fprintf('\n扫描 Level 2 指数数量:\n');
best_err = inf; best_M = 0;
for M_test = 2:min(20, ns)
    U_M = U(:,1:M_test); V_M = V(:,1:M_test); S_M = S(1:M_test,1:M_test);
    Z_M = S_M \ (U_M' * Y2mat * V_M);
    w_M = eig(Z_M);
    [~,idx] = sort(abs(w_M),'descend');
    w_M = w_M(idx);
    
    Y3m = zeros(N, M_test);
    for ii=1:M_test, Y3m(:,ii) = w_M(ii).^(0:N-1)'; end
    a_M = Y3m \ y2;
    at_M = log(w_M)/dt2;
    
    al2k_M = at_M.*T02./((1+1j*T02)*k);
    a2k_M = a_M.*exp(k*al2k_M);
    
    aa = [a1k(:);a2k_M(:)]; aal = [al1k(:);al2k_M(:)];
    G_M = compute_G(aa, aal, k, rho_test);
    err_M = abs(G_M-G_ref)/abs(G_ref)*100;
    
    fprintf('M=%2d: max|a_kz|=%.2e, abs(G)=%.2f, err=%6.2f%%\n', ...
        M_test, max(abs(a2k_M)), abs(G_M), err_M);
    
    if err_M < best_err
        best_err = err_M;
        best_M = M_test;
        best_aa = aa;
        best_aal = aal;
        best_G = G_M;
    end
end

fprintf('\n最优: M=%d, 误差=%.2f%%\n', best_M, best_err);
fprintf('G_dcim = %.4e %+.4ei (abs=%.4f)\n', real(best_G), imag(best_G), abs(best_G));
fprintf('G_ref  = %.4e %+.4ei (abs=%.4f)\n', real(G_ref), imag(G_ref), abs(G_ref));

function G = compute_G(aa, aal, k, rho)
    G = 0;
    for ii=1:length(aa)
        Rc = sqrt(rho^2-aal(ii)^2);
        G = G + aa(ii)*exp(-1j*k*Rc)/Rc;
    end
    G = G*1j/(2*pi);
end
