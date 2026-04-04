%% 最终方案：k=k_min + 转换后系数稳定性过滤
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

% 基准 - 所有分量
smgf2 = SpectralMGF();
smgf2.Initialize(lm, freq, [1,1,1,1,1], true, false);
smgf2.SetLayers(i_src, m_obs);
smgf2.SetSourcePoint(z_src); smgf2.SetObservationPoint(z_obs);
a_sw = 1.2*abs(lm.k_max);
qmgf = QuasistaticMGF2();
qmgf.Initialize(lm, freq, [1,1,1,1,1], false);
qmgf.SetLayers(i_src, m_obs);
K_q = qmgf.ComputeQMGF_Spatial(z_obs, z_src, rho_test);
orders = [0,1,1,0,0];
G_ref = zeros(1,5);
for cc=1:5
    G_int = SommerfeldIntegrator2.IntegrateSpectralNearField(smgf2, rho_test, cc, orders(cc), a_sw, false) ...
          + SommerfeldIntegrator2.IntegrateSpectralFarField(smgf2, rho_test, cc, orders(cc), a_sw, false);
    G_ref(cc) = G_int + K_q(cc);
end

%% DCIM 采样
N1 = 100; N2 = 100;
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

% 采样所有分量
K1_all = zeros(5,N1); K2_all = zeros(5,N2);
for ii=1:N1
    smgf.SetRadialWaveNumber(krho1(ii)); smgf.ComputeSpectralMGF();
    for cc=1:5, K1_all(cc,ii) = smgf.K(cc)*kz1(ii); end
end
for ii=1:N2
    smgf.SetRadialWaveNumber(krho2(ii)); smgf.ComputeSpectralMGF();
    for cc=1:5, K2_all(cc,ii) = smgf.K(cc)*kz2(ii); end
end

%% 对每个分量进行 DCIM + 扫描 M 选择最稳定的
fprintf('=== 每个分量独立优化 M ===\n');
for comp = 1:5
    K1 = K1_all(comp,:);
    K2 = K2_all(comp,:);
    
    % Level 1
    [at1, a1, ~] = GPOF(K1(:), t1(:), 50, 1e-4, 1e-16);
    a1k = a1.*exp(-T02*at1); al1k = at1./(1j*k);
    
    % Subtract L1
    K2m = K2;
    for ii=1:N2
        for jj=1:length(a1k)
            K2m(ii) = K2m(ii)-a1k(jj)*exp(-al1k(jj)*kz2(ii));
        end
    end
    
    % Level 2: SVD 准备
    y2 = K2m(:); t2v = t2(:); Nv = length(y2); Lv = 50;
    dt2 = t2v(2)-t2v(1);
    mr = Nv-Lv; nc = Lv;
    Y1mat = zeros(mr,nc); Y2mat = zeros(mr,nc);
    for ii=1:nc
        for jj=1:mr
            Y1mat(jj,ii) = y2(jj+ii-1);
            Y2mat(jj,ii) = y2(jj+ii);
        end
    end
    [U,S,V] = svd(Y1mat);
    sv = diag(S); ns = min(mr,nc);
    
    % 扫描 M，选择使 max|a_kz| 最小同时残差合理的 M
    best_err = inf; best_M = 2;
    for M_test = 2:min(15, ns)
        U_M = U(:,1:M_test); V_M = V(:,1:M_test); S_M = S(1:M_test,1:M_test);
        Z_M = S_M \ (U_M' * Y2mat * V_M);
        w_M = eig(Z_M);
        [~,idx] = sort(abs(w_M),'descend'); w_M = w_M(idx);
        
        Y3m = zeros(Nv,M_test);
        for ii=1:M_test, Y3m(:,ii) = w_M(ii).^(0:Nv-1)'; end
        a_M = Y3m \ y2;
        at_M = log(w_M)/dt2;
        
        al2k_M = at_M.*T02./((1+1j*T02)*k);
        a2k_M = a_M.*exp(k*al2k_M);
        
        aa = [a1k(:);a2k_M(:)]; aal = [al1k(:);al2k_M(:)];
        G_M = compute_G_comp(aa, aal, k, rho_test, comp);
        err_M = abs(G_M-G_ref(comp))/abs(G_ref(comp))*100;
        
        % 选择准则：最小误差（对于有基准的情况）
        % 实际使用时没有基准，需要用 max|a_kz| 作为代理
        if err_M < best_err
            best_err = err_M;
            best_M = M_test;
        end
    end
    
    % 也计算用 max|a_kz| 策略选择的结果
    best_err_stable = inf; best_M_stable = 2;
    for M_test = 2:min(15, ns)
        U_M = U(:,1:M_test); V_M = V(:,1:M_test); S_M = S(1:M_test,1:M_test);
        Z_M = S_M \ (U_M' * Y2mat * V_M);
        w_M = eig(Z_M);
        [~,idx] = sort(abs(w_M),'descend'); w_M = w_M(idx);
        
        Y3m = zeros(Nv,M_test);
        for ii=1:M_test, Y3m(:,ii) = w_M(ii).^(0:Nv-1)'; end
        a_M = Y3m \ y2;
        at_M = log(w_M)/dt2;
        
        al2k_M = at_M.*T02./((1+1j*T02)*k);
        a2k_M = a_M.*exp(k*al2k_M);
        res_M = norm(y2 - Y3m*a_M)/norm(y2);
        
        % 策略：选择 max|a_kz| < 10 且残差最小的
        if max(abs(a2k_M)) < 10 && res_M < best_err_stable
            best_err_stable = res_M;
            best_M_stable = M_test;
        end
    end
    
    % 输出两种策略
    % 最优M（oracle）
    U_M = U(:,1:best_M); V_M = V(:,1:best_M); S_M = S(1:best_M,1:best_M);
    Z_M = S_M \ (U_M' * Y2mat * V_M);
    w_M = eig(Z_M); [~,idx] = sort(abs(w_M),'descend'); w_M = w_M(idx);
    Y3m = zeros(Nv,best_M); for ii=1:best_M, Y3m(:,ii) = w_M(ii).^(0:Nv-1)'; end
    a_M = Y3m \ y2; at_M = log(w_M)/dt2;
    al2k_M = at_M.*T02./((1+1j*T02)*k); a2k_M = a_M.*exp(k*al2k_M);
    aa = [a1k(:);a2k_M(:)]; aal = [al1k(:);al2k_M(:)];
    G_oracle = compute_G_comp(aa, aal, k, rho_test, comp);
    
    % 稳定M
    U_M = U(:,1:best_M_stable); V_M = V(:,1:best_M_stable); S_M = S(1:best_M_stable,1:best_M_stable);
    Z_M = S_M \ (U_M' * Y2mat * V_M);
    w_M = eig(Z_M); [~,idx] = sort(abs(w_M),'descend'); w_M = w_M(idx);
    Y3m = zeros(Nv,best_M_stable); for ii=1:best_M_stable, Y3m(:,ii) = w_M(ii).^(0:Nv-1)'; end
    a_M = Y3m \ y2; at_M = log(w_M)/dt2;
    al2k_M = at_M.*T02./((1+1j*T02)*k); a2k_M = a_M.*exp(k*al2k_M);
    aa_s = [a1k(:);a2k_M(:)]; aal_s = [al1k(:);al2k_M(:)];
    G_stable = compute_G_comp(aa_s, aal_s, k, rho_test, comp);
    
    err_oracle = abs(G_oracle-G_ref(comp))/abs(G_ref(comp))*100;
    err_stable = abs(G_stable-G_ref(comp))/abs(G_ref(comp))*100;
    
    fprintf('comp%d: oracle M=%2d err=%6.2f%% | stable(|a|<10) M=%2d err=%6.2f%%\n', ...
        comp, best_M, err_oracle, best_M_stable, err_stable);
end

function G = compute_G_comp(aa, aal, k, rho, comp)
    G = 0;
    rho_sq = rho^2;
    for ii=1:length(aa)
        Rc = sqrt(rho_sq-aal(ii)^2);
        if comp==1||comp==4||comp==5
            G = G + aa(ii)*exp(-1j*k*Rc)/Rc;
        else
            G = G + rho*aa(ii)*(1+1j*k*Rc)*exp(-1j*k*Rc)/(Rc^3);
        end
    end
    G = G*1j/(2*pi);
end
