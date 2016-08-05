function [lambda] = GeometricMean_MatlabMethod(fxy,m,n_k)
% Compute the geometric mean of the entries of the Cauchy like matrix
% C_{n-k}(f(x,y))
%
% % Inputs.
%
% fxy : Coefficients of polynomial f(x,y)
%
% m : Total degree of f(x,y)
%
% n_k : Total degree of v_{k}(x,y)

% Build the diagonal matrix D^{-1}
D = BuildD(m,n_k);

% Build the matrix T_{n-k}(f(x,y))
T1 = BuildT1(fxy,m,n_k);

% Build the matrix Q
Q1 = BuildQ1(n_k);

% Build the matrix D_{m+n-k}^{-1}T_{n-k}(f(x,y))*Q_{n-k}
DT1Q1 = D * T1 * Q1;

% Compute geometric mean
lambda = geomean(abs(DT1Q1(DT1Q1~=0)));


end