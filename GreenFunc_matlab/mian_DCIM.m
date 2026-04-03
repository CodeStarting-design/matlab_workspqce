%% main.m - Test Driver for Multilayer Green's Function
clear; clc;

% 1. 设置仿真参数
freq = 30e9; % 30 GHz
omega = 2 * pi * freq;
lambda0 = 2*pi/(omega*sqrt(Constants.eps0*Constants.mu0));

% 2. 定义层结构 (对应 layers.txt)
% 注意单位转换: layers.txt 中单位是 mm
unit = 1e-3;

lm = LayerManager();

% 添加介质层 (从上到下)
% L1: zmin=1.1, h=0.7 -> zmax=1.8
lm.AddLayer(1.1*unit, 1.8*unit, 2.1, 1, 0); 

% L2: zmin=0.8, h=0.3 -> zmax=1.1
lm.AddLayer(0.8*unit, 1.1*unit, 12.5, 1, 0);

% L3: zmin=0.3, h=0.5 -> zmax=0.8
lm.AddLayer(0.3*unit, 0.8*unit, 9.8, 1, 0);

% L4: zmin=0, h=0.3 -> zmax=0.3
lm.AddLayer(0.0*unit, 0.3*unit, 8.6, 1, 0);


% 设置半无限空间
% Top: Air, Bot: PEC (sigma = -1 in config -> isPEC=true)
lm.SetHalfspaces(1.0, 1.0, 0.0, ...   % Top: eps, mu, sigma
                 1.0, 1.0, 0.0, ...   % Bot: eps, mu, sigma (PEC ignored here)
                 false, true);        % isPEC_Top, isPEC_Bot

% 处理层参数
lm.ProcessLayers(freq);

% 打印层信息以确认
disp('Layers processed.');

% 设置DCIM参数
components = [1,1,1,1,1];  % Gxx和Gphi

% 创建DCIM对象
dcim = DCIM_Integrated();
dcim.extract_quasistatic = true;  % 启用准静态项提取
dcim.N1 = 100;  % 增加采样点数
dcim.N2 = 100;

% 初始化
dcim.Initialize(lm, freq, components);

% 4. 设置源和观测点
x_src = 0.0;
y_src = 0.0;
z_src = 0.4 * unit; % 0.4 mm
% x_obs = 0.0;
y_obs = 0.0;
z_obs = 1.4 * unit; % 1.4 mm (Note: This is in L1? Check layers)
% L4: 0-0.3, L3: 0.3-0.8, L2: 0.8-1.1, L1: 1.1-1.8. 
% z_src=0.4 is inside L3. z_obs=1.4 is inside L1.

i_src = lm.FindLayer(z_src);
m_obs = lm.FindLayer(z_obs);

% 生成DCIM复镜像
fprintf('\n生成DCIM复镜像...\n');
dcim.GenerateImages(i_src, m_obs, z_src, z_obs);

% 5. 扫描位置并计算
x_min = 1.6e-4 * lambda0;
x_max = 1.6e1 * lambda0;
Nx = 500;
x_vec = logspace(log10(x_min), log10(x_max), Nx);
G_results = zeros(Nx, 10); % 9 dyadic + 1 scalar

%% 计算DCIM结果
fprintf('Starting computation loop...\n');
for ii = 1:Nx
    x_obs = x_vec(ii);
    x_diff = x_obs - x_src;
    y_diff = y_obs - y_src;
    [G_results(ii,1:9), G_results(ii,10)] = dcim.ComputeSpatialMGF(x_diff, y_diff, z_obs, z_src);
end
fprintf('Computation done.\n');

% 6. 绘图 (比如 Gxx 和 Gphi)
figure;
test = readtable('MGFdata.txt');
loglog(table2array(test(:,1))/lambda0, table2array(test(:,2)), 'b-','LineWidth', 1.5);
hold on
loglog(table2array(test(:,1))/lambda0, table2array(test(:,4)), 'b-','LineWidth', 1.5);
hold on
loglog(table2array(test(:,1))/lambda0, table2array(test(:,8)), 'b-','LineWidth', 1.5);
hold on
loglog(table2array(test(:,1))/lambda0, table2array(test(:,10)), 'b-','LineWidth', 1.5);
hold on
loglog(table2array(test(:,1))/lambda0, table2array(test(:,11)), 'b-','LineWidth', 1.5);
hold on

loglog(x_vec/lambda0, abs(G_results(:, 1)), 'r--', 'LineWidth', 1.5); hold on;
loglog(x_vec/lambda0, abs(G_results(:, 3)), 'r--', 'LineWidth', 1.5); hold on;
loglog(x_vec/lambda0, abs(G_results(:, 7)), 'r--', 'LineWidth', 1.5); hold on;
loglog(x_vec/lambda0, abs(G_results(:, 9)), 'r--', 'LineWidth', 1.5); hold on;
loglog(x_vec/lambda0, abs(G_results(:, 10)), 'r--', 'LineWidth', 1.5);
xlabel('rho / lambda_0');
ylabel('Magnitude');
legend('Gxx','Gxz','Gzx','Gzz','Gphi');
grid on;
title('MGF Components vs Distance');
