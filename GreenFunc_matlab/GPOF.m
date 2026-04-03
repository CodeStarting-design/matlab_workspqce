function [exponents, coefficients, residual] = GPOF(y, t, L, tol_svd, tol_eig)
%GPOF Generalized Pencil of Function for complex-valued data
%
% Fits a function y(t) as a sum of complex exponentials:
%   y(t) = sum_{i=1}^{M} a_i * exp(alpha_i * t)
%
% This implementation handles complex-valued functions, following
% the strata C++ implementation.
%
% Inputs:
%   y       - Sampled function values (column vector, complex-valued)
%   t       - Sample points (column vector, real-valued, uniformly spaced)
%   L       - Pencil parameter (default: N/2)
%   tol_svd - Tolerance for singular value selection (default: 1e-4)
%   tol_eig - Tolerance for eigenvalue selection (default: 1e-16)
%
% Outputs:
%   exponents    - Complex exponents alpha_i
%   coefficients - Complex coefficients a_i
%   residual     - Fitting residual
%
% Reference:
%   strata/src/DCIM.cpp RunGPOF()
%   Y. Hua and T. K. Sarkar, "Matrix pencil method for estimating
%   parameters of exponentially damped/undamped sinusoids in noise",
%   IEEE Trans. Acoust., Speech, Signal Process., vol. 38, no. 5, 1990.

if nargin < 3 || isempty(L)
    L = [];
end

if nargin < 4 || isempty(tol_svd)
    tol_svd = 1e-4;
end

if nargin < 5 || isempty(tol_eig)
    tol_eig = 1e-16;
end

% 确保y和t是列向量
y = y(:);
t = t(:);

N = length(y);
if N ~= length(t)
    error('y and t must have the same length');
end

if isempty(L)
    L = floor(N / 2);
end

if L >= N
    error('Pencil parameter L must be less than length(y)');
end

% 检查采样是否均匀
dt = t(2) - t(1);
if any(abs(diff(t) - dt) > 1e-10 * abs(dt))
    warning('GPOF: Sampling is not uniform. Results may be inaccurate.');
end

% Construct Y1 and Y2 matrices (strata style, column-major)
m = N - L;  % rows
n = L;      % columns

Y1 = zeros(m, n);
Y2 = zeros(m, n);

for ii = 1:n        % columns
    for jj = 1:m    % rows
        Y1(jj, ii) = y(jj + ii - 1);
        Y2(jj, ii) = y(jj + ii);
    end
end

% SVD of Y1
[U, S, V] = svd(Y1, 'econ');
singular_values = diag(S);

% Determine effective rank using tolerance (strata style)
M = length(singular_values);
for M = 1:length(singular_values)
    if M < length(singular_values) && singular_values(M+1) < tol_svd * singular_values(1)
        break;
    end
end

% Truncate to rank M
U_M = U(:, 1:M);
V_M = V(:, 1:M);
S_M = S(1:M, 1:M);

% Compute Z matrix: Z = S_M^{-1} * U_M' * Y2 * V_M
Z = S_M \ (U_M' * Y2 * V_M);

% Compute eigenvalues of Z
w = eig(Z);

% Filter eigenvalues by magnitude
[~, idx] = sort(abs(w), 'descend');
w = w(idx);

nwt = 1;
for nwt = 1:length(w)
    if nwt < length(w) && abs(w(nwt+1)) < tol_eig * abs(w(1))
        break;
    end
end


w = w(1:nwt);

% Note: strata does NOT filter eigenvalues by |w| < 1
% It uses all eigenvalues that pass the tol_eig test
% We follow strata's approach

% Set up the system of equations for coefficients
Y3 = zeros(N, nwt);
for ii = 1:nwt
    Y3(:, ii) = w(ii).^(0:N-1)';
end

% Solve least squares: Y3 * a = y
a = Y3 \ y;

% Extract exponents and coefficients
coefficients = a;
exponents = log(w) / dt;

% Compute residual
y_fit = zeros(N, 1);
for ii = 1:nwt
    y_fit = y_fit + coefficients(ii) * exp(exponents(ii) * t);
end

residual = norm(y - y_fit) / norm(y);

end
