classdef Constants
    properties (Constant)
        % 基本数学与物理常数
        j = 1i; % MATLAB 使用 1i
        pi = pi;
        
        eps0 = 8.854187817e-12;       % 真空介电常数
        mu0 = 4.0 * pi * 1.0e-7;      % 真空磁导率
        c0 = 1.0 / sqrt(8.854187817e-12 * 4.0 * pi * 1.0e-7); % 光速
        eta0 = sqrt(4.0 * pi * 1.0e-7 / 8.854187817e-12);   % 波阻抗
        
        EulerGamma = 0.5772156649015328606065; % 欧拉常数
        inf = Inf;
    end
end