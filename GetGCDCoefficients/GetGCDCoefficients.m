function dxy = GetGCDCoefficients(fxy,gxy,uxy,vxy,m,n,t)
% Compute coefficients of polynomial d(x,y)

% % Inputs.
%
% fxy : Coefficients of polynomial f(x,y)
% 
% gxy : Coefficients of polynomial g(x,y)
%
% uxy : Coefficients of polynomial u(x,y)
%
% vxy : Coefficients of polynomial v(x,y)
%
% m : Total degree of f(x,y)
%
% n : Total degree of g(x,y)
%
% t : Total degree of d(x,y)

% % 
% Build 

% Build matrix H
H1 = BuildD(m-t,t);
H2 = BuildD(n-t,t);
H = blkdiag(H1,H2);

% Build matrix C_{t}(u)
C1 = BuildT1(uxy,m-t,t);

% Build matrix C_{t}(v)
C2 = BuildT1(vxy,n-t,t);

C = [C1 ; C2];

% Build Matrix of binomial corresponding to d(x,y)
Q = BuildQ1(t);

% Build RHS vector
v_fxy = GetAsVector(fxy);
nCoefficients_fxy = nchoosek(m+2,2);
v_fxy = v_fxy(1:nCoefficients_fxy);

v_gxy = GetAsVector(gxy);
nCoefficients_gxy = nchoosek(n+2,2);
v_gxy = v_gxy(1:nCoefficients_gxy);

rhs_vec = [v_fxy;v_gxy];

% Build System
x_ls = SolveAx_b(H*C*Q,rhs_vec);

try
nZeros_dxy = nchoosek(t+1,2);
catch
    nZeros_dxy = 0;
end

v_dxy = [...
    x_ls;
    zeros(nZeros_dxy,1)
    ];

dxy = GetAsMatrix(v_dxy,t,t);

end