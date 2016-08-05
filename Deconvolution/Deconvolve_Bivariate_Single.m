function hxy = Deconvolve_Bivariate_Single(fxy,gxy)
%
% Deconvolve the polynomials f(x,y) and g(x,y) to obtain h(x,y)
%
% % Inputs.
%
% fxy : Coefficients of polynomial f(x,y)
%
% gxy : Coefficients of polynomial g(x,y)
%
% m : Degree of polynomial f(x,y)
%
% n : Degree of polynomial g(x,y)

[m1,m2] = GetDegree(fxy);
[n1,n2] = GetDegree(gxy);

m = m1;
n = n1;

% Build the matrix D^{-1}
D = BuildD(n,m-n);

% Build the matrix T_{m-n}(g)
T1 = BuildT1(gxy,n,m-n);

% Build the matrix Q_{m-n}
Q = BuildQ1(m-n);

% Build the matrix D^{-1}_{m}*T_{m-n}(g)Q_{m-n}
DTQ = D*T1*Q;

% Build the rhs coefficient vector
f = GetAsVector(fxy);
nCoefficients_fxy = nchoosek(m+2,2);
v_fxy = f(1:nCoefficients_fxy);

% Get coefficients of h(x,y)
x_ls = SolveAx_b(DTQ,v_fxy);

% Get number of zeros in the coefficient matrix for polynomial h(x,y)
if (m-n+1 > 1)
    nZeros_hxy = nchoosek(m-n+1,2);
else
    nZeros_hxy = 0;
end
v_hxy = ...
    [
    x_ls;
    zeros(nZeros_hxy,1);
    ];

hxy = GetAsMatrix(v_hxy,m-n,m-n);



end

