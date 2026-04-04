function [exponents, coefficients, residual] = GPOF(y, t, L, tol_svd, tol_eig)
%GPOF Generalized Pencil of Function for complex-valued data
%
% Fits a function y(t) as a sum of complex exponentials:
%   y(t) = sum_{i=1}^{M} a_i * exp(alpha_i * t)
%
% This implementation follows the strata C++ RunGPOF() exactly:
%   - Full SVD (not economy)
%   - Zero out D_inv entries instead of truncating matrix dimensions
%   - Eigenvalues of full ns x ns Z matrix
%
% Reference: strata/src/DCIM.cpp RunGPOF()

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

% Construct Y1 and Y2 matrices (strata style)
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

% Full SVD of Y1 (C++ uses 'A','A' = full SVD)
[U, S, V] = svd(Y1);
singular_values = diag(S);
ns = min(m, n);  % number of singular values

% Determine truncation point (strata style)
nst = ns;  % default: keep all
for ii = 2:ns
    if singular_values(ii) < tol_svd * singular_values(1)
        nst = ii - 1;
        break;
    end
end

% Build D_inv: ns x ns diagonal matrix, zero out small singular values
% (C++ style: don't truncate matrix, just zero out)
D_inv = zeros(ns, ns);
for ii = 1:nst
    D_inv(ii, ii) = 1.0 / singular_values(ii);
end
% entries nst+1:ns remain zero

% Compute Z = D_inv * U_ns^H * Y2 * V_ns  (all ns x ns)
% U_ns = first ns columns of U (m x ns)
% V_ns = first ns columns of V (n x ns)
U_ns = U(:, 1:ns);
V_ns = V(:, 1:ns);

c1 = U_ns' * Y2 * V_ns;   % ns x ns
Z = D_inv * c1;            % ns x ns

% Compute eigenvalues of Z (full ns x ns matrix, like C++)
w = eig(Z);

% Sort by magnitude (descending)
[~, idx] = sort(abs(w), 'descend');
w = w(idx);

% Filter by tol_eig (strata style)
nwt = ns;
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
