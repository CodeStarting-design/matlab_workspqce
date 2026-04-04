function [exponents, coefficients, residual] = GPOF(y, t, L, tol_svd, tol_eig)
%GPOF Generalized Pencil of Function for complex-valued data
%
% Fits y(t) = sum a_i * exp(alpha_i * t)
%
% Uses TLS-based Matrix Pencil Method following strata C++ RunGPOF(),
% with improved numerical stability via SVD-filtered generalized
% eigenvalue approach.

if nargin < 3 || isempty(L), L = []; end
if nargin < 4 || isempty(tol_svd), tol_svd = 1e-4; end
if nargin < 5 || isempty(tol_eig), tol_eig = 1e-16; end

y = y(:);
t = t(:);
N = length(y);
if N ~= length(t), error('y and t must have the same length'); end
if isempty(L), L = floor(N / 2); end
if L >= N, error('Pencil parameter L must be less than length(y)'); end

dt = t(2) - t(1);

% Construct Y1 and Y2 matrices
m = N - L;  % rows
n = L;      % columns

Y1 = zeros(m, n);
Y2 = zeros(m, n);
for ii = 1:n
    for jj = 1:m
        Y1(jj, ii) = y(jj + ii - 1);
        Y2(jj, ii) = y(jj + ii);
    end
end

% Full SVD of Y1
[U, S, V] = svd(Y1);
singular_values = diag(S);
ns = min(m, n);

% Determine truncation point
nst = ns;
for ii = 2:ns
    if singular_values(ii) < tol_svd * singular_values(1)
        nst = ii - 1;
        break;
    end
end

% Truncated SVD approach: Z = S_nst^{-1} * U_nst^H * Y2 * V_nst
% This gives an nst x nst matrix whose eigenvalues are the signal poles
U_nst = U(:, 1:nst);
V_nst = V(:, 1:nst);
S_nst_inv = diag(1.0 ./ singular_values(1:nst));

Z = S_nst_inv * (U_nst' * Y2 * V_nst);

% Compute eigenvalues
w = eig(Z);

% Sort by magnitude (descending)
[~, idx] = sort(abs(w), 'descend');
w = w(idx);

% Filter by tol_eig
nwt = length(w);
for ii = 2:length(w)
    if abs(w(ii)) < tol_eig * abs(w(1))
        nwt = ii - 1;
        break;
    end
end
w = w(1:nwt);

% Set up system for coefficients
Y3 = zeros(N, nwt);
for ii = 1:nwt
    Y3(:, ii) = w(ii).^(0:N-1)';
end

% Solve least squares
a = Y3 \ y;

% Extract results
coefficients = a;
exponents = log(w) / dt;

% Compute residual
y_fit = Y3 * a;
residual = norm(y - y_fit) / norm(y);

end
