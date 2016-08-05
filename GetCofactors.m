function [uxy,vxy] = GetCofactors(fxy,gxy,t)

[m1,m2] = GetDegree(fxy);
m = m1;

[n1,n2] = GetDegree(gxy);
n = n1;

% Build the Sylvester matrix
D = BuildD(m,n-t);
T1 = BuildT1(fxy,m,n-t);
T2 = BuildT1(gxy,n,m-t);
Q = BuildQ(m,n,t);
Sk = D* [T1 T2] * Q;

% Remove optimal column
opt_col = GetOptimalColumn(Sk);

% Get A_{k}, the matrix S_{k} with opt column removed
Ak = Sk;
Ak(:,opt_col) = [];

% Get c_{k}, the column removed from S_{k}
ck = Sk(:,opt_col);

% Solve Ax = b
x_ls = SolveAx_b(Ak,ck);

% Insert '-1' into the opt_col position of x_ls
x = ...
    [...
        x_ls(1:opt_col-1) 
        -1
        x_ls(opt_col:end)
    ];

% Split x to obtain v(x,y) and u(x,y)
nCoefficients_vxy = nchoosek(n-t+2,2);


% Get Coefficients of u(x,y) and v(x,y)
vxy = x(1:nCoefficients_vxy);
uxy = -1 .* x(nCoefficients_vxy+1:end);

% Include zero coefficients (so that u(x,y) and v(x,y) can be placed into
% (m-t+1) x (m-t+1) and (n-t+1) x (n-t+1) matrices respectively, where the
% coefficients take the upper left triangle.
try
nZeros_uxy = nchoosek(m-t+1,2);
catch
    nZeros_uxy = 0;
end
try
nZeros_vxy = nchoosek(n-t+1,2);
catch
    nZeros_vxy = 0;
end
uxy = [uxy ; zeros(nZeros_uxy,1)];
vxy = [vxy ; zeros(nZeros_vxy,1)];

uxy = GetAsMatrix(uxy,m-t,m-t);
vxy = GetAsMatrix(vxy,n-t,n-t);
    



end
