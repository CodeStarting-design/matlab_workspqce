%% 分析转换后系数稳定性
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

for tol_svd = [1e-4, 1e-6]
    fprintf('\n======== tol_svd = %.0e ========\n', tol_svd);
    
    [at1, a1, ~] = GPOF(K1(:), t1(:), 50, tol_svd, 1e-16);
    a1k = a1.*exp(-T02*at1); al1k = at1./(1j*k);
    
    fprintf('L1: %d exp\n', length(at1));
    for jj=1:length(a1k)
        fprintf('  a_t=%.3e, alpha_t=%.3e%+.3ei -> a_kz=%.3e, alpha_kz=%.3e%+.3ei\n', ...
            abs(a1(jj)), real(at1(jj)), imag(at1(jj)), abs(a1k(jj)), real(al1k(jj)), imag(al1k(jj)));
    end
    
    K2m = K2;
    for ii=1:N2
        for jj=1:length(a1k)
            K2m(ii) = K2m(ii)-a1k(jj)*exp(-al1k(jj)*kz2(ii));
        end
    end
    
    [at2, a2, res2] = GPOF(K2m(:), t2(:), 50, tol_svd, 1e-16);
    al2k = at2.*T02./((1+1j*T02)*k);
    a2k = a2.*exp(k*al2k);
    
    fprintf('L2: %d exp, res=%.4e\n', length(at2), res2);
    for jj=1:length(a2k)
        fprintf('  a_t=%.3e, alpha_t=%.3e%+.3ei -> a_kz=%.3e, alpha_kz=%.3e%+.3ei\n', ...
            abs(a2(jj)), real(at2(jj)), imag(at2(jj)), abs(a2k(jj)), real(al2k(jj)), imag(al2k(jj)));
    end
    
    % Sommerfeld identity贡献分析
    aa = [a1k(:);a2k(:)]; aal = [al1k(:);al2k(:)];
    fprintf('\nSommerfeld identity (rho=%.2e):\n', rho_test);
    G_total = 0;
    for ii=1:length(aa)
        Rc = sqrt(rho_test^2-aal(ii)^2);
        contrib = aa(ii)*exp(-1j*k*Rc)/Rc * 1j/(2*pi);
        G_total = G_total + contrib;
        fprintf('  img%2d: |a|=%.3e, alpha=%.3e%+.3ei, |Rc|=%.3e, |contrib|=%.3e\n', ...
            ii, abs(aa(ii)), real(aal(ii)), imag(aal(ii)), abs(Rc), abs(contrib));
    end
    fprintf('G_dcim = %.4e %+.4ei (abs=%.4f)\n', real(G_total), imag(G_total), abs(G_total));
    fprintf('G_ref  = %.4e %+.4ei (abs=%.4f)\n', real(G_ref), imag(G_ref), abs(G_ref));
    fprintf('误差 = %.2f%%\n', abs(G_total-G_ref)/abs(G_ref)*100);
end
