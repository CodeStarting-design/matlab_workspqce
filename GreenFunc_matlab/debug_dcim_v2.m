%% 测试新 GPOF 实现
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
krho1 = sqrt(k^2-kz1.^2);
krho1 = abs(real(krho1))+1j*abs(imag(krho1));

t2 = linspace(1e-2, T02, N2);
kz2 = k*(-1j*t2+(1-t2/T02));
krho2 = sqrt(k^2-kz2.^2);
krho2 = abs(real(krho2))+1j*abs(imag(krho2));

K1 = zeros(1,N1);
for ii=1:N1
    smgf.SetRadialWaveNumber(krho1(ii));
    smgf.ComputeSpectralMGF();
    K1(ii) = smgf.K(comp)*kz1(ii);
end

[at1, a1, res1] = GPOF(K1(:), t1(:), 50, 1e-4, 1e-16);
fprintf('L1: %d exp, res=%.4e\n', length(at1), res1);
a1k = a1.*exp(-T02*at1); al1k = at1./(1j*k);

K2 = zeros(1,N2);
for ii=1:N2
    smgf.SetRadialWaveNumber(krho2(ii));
    smgf.ComputeSpectralMGF();
    K2(ii) = smgf.K(comp)*kz2(ii);
end

K2m = K2;
for ii=1:N2
    for jj=1:length(a1k)
        K2m(ii) = K2m(ii)-a1k(jj)*exp(-al1k(jj)*kz2(ii));
    end
end

[at2, a2, res2] = GPOF(K2m(:), t2(:), 50, 1e-4, 1e-16);
fprintf('L2: %d exp, res=%.4e\n', length(at2), res2);
al2k = at2.*T02./((1+1j*T02)*k);
a2k = a2.*exp(k*al2k);

aa = [a1k(:);a2k(:)]; aal = [al1k(:);al2k(:)];
fprintf('总镜像: %d\n', length(aa));

G_dcim = 0;
for ii=1:length(aa)
    Rc = sqrt(rho_test^2-aal(ii)^2);
    G_dcim = G_dcim + aa(ii)*exp(-1j*k*Rc)/Rc;
end
G_dcim = G_dcim*1j/(2*pi);

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

fprintf('\nG_dcim = %.6e %+.6ei (abs=%.4f)\n', real(G_dcim), imag(G_dcim), abs(G_dcim));
fprintf('G_ref  = %.6e %+.6ei (abs=%.4f)\n', real(G_ref), imag(G_ref), abs(G_ref));
fprintf('误差 = %.2f%%\n', abs(G_dcim-G_ref)/abs(G_ref)*100);

% 实轴验证
fprintf('\n实轴验证:\n');
krho_tests = [10, 500, 800, 1000, 1500, 2000];
for kk = 1:length(krho_tests)
    kr = krho_tests(kk)+0.001j;
    kz_t = sqrt(k^2-kr^2);
    if imag(kz_t)>0, kz_t = complex(real(kz_t),-abs(imag(kz_t))); end
    if real(kz_t)<0, kz_t = complex(abs(real(kz_t)),imag(kz_t)); end
    smgf.SetRadialWaveNumber(kr);
    smgf.ComputeSpectralMGF();
    K_dir = smgf.K(comp)*kz_t;
    K_dcim_v = 0;
    for jj=1:length(aa)
        K_dcim_v = K_dcim_v+aa(jj)*exp(-aal(jj)*kz_t);
    end
    fprintf('  krho=%5.0f: err=%.1f%%\n', real(kr), abs(K_dcim_v-K_dir)/abs(K_dir)*100);
end
