function [fxy_lr,gxy_lr,uxy_lr,vxy_lr] = STLN(fxy,gxy,k)
% Compute the low rank approximation of the Sylvester matrix S_{k}(f,g)
% by method of Structured Total Least Norm STLN.
%
% % Inputs.
%
% fxy : Coefficients of polynomial f(x,y)
%
% gxy : Coefficients of polynomial g(x,y)
%
% k : Degree of common divisor d(x,y) and index of the Sylvester matrix
%     S_{k} whose low rank approximation is considered.
%
% % Outputs.
%
% fxy_lr : Coefficients f(x,y) obtained from the low rank approx of S_{k}(f,g)
%
% gxy_lr : Coefficients g(x,y) obtained from the low rank approx of S_{k}(f,g)
%
% uxy_lr : Coefficients u(x,y) obtained from the low rank approx of S_{k}(f,g)
%
% vxy_lr : Coefficients v(x,y) obtained from the low rank approx of S_{k}(f,g)


% Get the degree of f(x,y)
[m,~] = GetDegree(fxy);

% Get the degree of g(x,y)
[n,~] = GetDegree(gxy);

% Get number of coefficients in f(x,y)
nCoeffs_fxy = nchoosek(m+2,2);

% Get the number of zeros in the matrix containing f(x,y)
nZeros_fxy = nchoosek(m+1,2);

% Get number of coefficients in g(x,y)
nCoeffs_gxy = nchoosek(n+2,2);

% Get the number of zeros in the matrix containing g(x,y)
nZeros_gxy = nchoosek(n+1,2);

% Get number of coefficients in u(x,y)
nCoeffs_uxy = nchoosek(m-k+2,2);

% Get the number of zeros in the matrix of u(x,y)
nZeros_uxy = nchoosek(m-k+1,2);

% Get number of coefficients in v(x,y)
nCoeffs_vxy = nchoosek(n-k+2,2);

% Get number of zeros in the matrix of v(x,y)
nZeros_vxy = nchoosek(n-k+1,2);

% Build the kth Sylvester Subresultant matrix S_{k}(f,g)
% Build the diagonal matrix D^{-1}
D = BuildD(m,n-k);

% Build the matrix T_{n-k}(f)
T1_fx = BuildT1(fxy,m,n-k);

% Build the matrix T_{m-k}(g)
T1_gx = BuildT1(gxy,n,m-k);

% Build the matrix Q_{k}
Q = BuildQ(m,n,k);

% Build the kth Sylvester subresultant matrix.
Sk_fg = D*[T1_fx T1_gx] * Q;

% Get index of optimal column for removal from S_{k}(f,g)
idx_col = GetOptimalColumn(Sk_fg);

% Remove optimal column from S_{k}(f,g)
Ak_fg = Sk_fg;
Ak_fg(:,idx_col) = [];

% Get c_{k} optimal column removed from S_{k}(f,g)
ck = Sk_fg(:,idx_col);

% Get least squares solution x_ls
x_ls = SolveAx_b(Ak_fg,ck);

% Get the residual vector associated with LS solution
res_vec = (ck) - (Ak_fg*x_ls);

% Build P_{t} - The matrix P_{t} such that P_{t}z = ht or P_{t}[f;g] = c_{t}
DPQ = BuildP_STLN(m,n,k,idx_col);

% Obtain the vector \hat{x} which contains x_ls with a zero inserted into
% the position of the optimal column.
xk = [x_ls(1:idx_col-1) ; 0 ; x_ls(idx_col:end)];

% Build matrix Y_{k}
DYQ = BuildY_STLN(xk,m,n,k);

% Get matrix of coefficients of polynomial z_{f}(x,y)
zf = zeros(nCoeffs_fxy,1);
vec_zf = [zf ; zeros(nZeros_fxy,1)];
mat_zf = GetAsMatrix(vec_zf,m,m);

% Get matrix of coefficients of polynomial z_{g}(x,y)
zg = zeros(nCoeffs_gxy,1);
vec_zg = [zg ; zeros(nZeros_gxy,1)];
mat_zg = GetAsMatrix(vec_zg,n,n);

% Get vector z = [z_{f} z_{g}]
vec_z = [zf;zg];

% Build the matrix T_{n-k}(z_{f})
T1_zf = BuildT1(mat_zf,m,n-k);

% Build the matrix T_{m-k}(z_{g})
T1_zg = BuildT1(mat_zg,n,m-k);

% Build the matrix S_{k}(z_{f},z_{g})
Sk_zfzg = D*[T1_zf T1_zg] * Q;
Ak_zfzg = Sk_zfzg;
Ak_zfzg(:,idx_col) = [];

% Build matrices H_{z} and H_{x}
H_z = DYQ - DPQ;
H_x = Ak_fg + Ak_zfzg;

% Build matrix C.
C = [H_z H_x];


% Build the matrix E for LSE Problem
nEntries_x = nchoosek(m-k+2,2) + nchoosek(n-k+2,2) - 1;
nEntries_z = nchoosek(m+2,2) + nchoosek(n+2,2);

E = eye(nEntries_x + nEntries_z);
%E = blkdiag(eye(nEntries_x), zeros(nEntries_z));

% Define the starting vector for the iterations for the LSE problem.
start_point     =   ...
    [...
    vec_z;...
    x_ls;
    ];

% Set yy to be the vector which stores all cummulative perturbations.
yy = start_point;

% Set the initial value of vector f to be zero
f = -(yy - start_point);

% Initialise the iteration counter
ite = 1;

% Set the termination criterion
condition(ite) = norm(res_vec)/norm(ck) ;


global SETTINGS
while condition(ite) > SETTINGS.STLN_MAX_ERROR && ite < SETTINGS.STLN_MAX_ITERATIONS
    
    % Increment interation counter
    ite = ite + 1;
    
    % Get small petrubations by LSE
    y_lse = LSE(E,f,C,res_vec);
    
    % Increment cummulative peturbations
    yy = yy + y_lse;
    
    % update z and x_ls
    delta_z = y_lse(1:nEntries_z);
    delta_xls = y_lse(nEntries_z + 1:end);
    
    % Get z1 and z2
    delta_zf = delta_z(1:nCoeffs_fxy);
    delta_zg = delta_z(nCoeffs_fxy+1:nCoeffs_fxy+nCoeffs_gxy);
    
    % Update vectors z, z_{f} and z_{g}
    vec_z = vec_z + delta_z;
    zf =  zf + delta_zf;
    zg =  zg + delta_zg;
    
    % % Build matrix S_{t}(z1,z2)
    vec_zf = [zf ; zeros(nZeros_fxy,1)];
    mat_zf = GetAsMatrix(vec_zf,m,m);
    
    vec_zg = [zg ; zeros(nZeros_gxy,1)];
    mat_zg = GetAsMatrix(vec_zg,n,n);
    
    T1_zf = BuildT1(mat_zf,m,n-k);
    T1_zg = BuildT1(mat_zg,n,m-k);
    
    % Build Sylvester matrix S_{t}(zf,zg)
    Sk_zfzg = D*[T1_zf T1_zg] * Q;
    
    % Get the matrix A_{t}(zf,zg)
    Ak_zfzg = Sk_zfzg;
    Ak_zfzg(:,idx_col) = [];
    
    % Get the vector h_{t}
    hk = Sk_zfzg(:,idx_col);
    
    % Get x_{1} and x_{2}
    
    % Obtain the vector \hat{x} which contains x_ls with a zero inserted into
    % the position of the optimal column.
    x_ls = x_ls + delta_xls;
    
    
    xk = [x_ls(1:idx_col-1) ; 0 ; x_ls(idx_col:end)];
    
    
    % Build the matrix Y_{k}
    DYQ = BuildY_STLN(xk,m,n,k);
    
    % % % Build the matrix C in the LSE Problem.
    
    % Build H_z
    H_z = DYQ - DPQ;
    
    % Build Hx
    H_x = (Ak_fg + Ak_zfzg);
    
    % Build C
    C = [H_z H_x];
    
    % Get the residual vector
    res_vec = (ck + hk) - ((Ak_fg + Ak_zfzg)*x_ls);
    
    
    % Update f - used in LSE Problem.
    f = -(yy - start_point);
    
    % Update the termination criterion
    condition(ite) = norm(res_vec)./ norm (ck + hk) ;
    
end

LineBreakLarge()
fprintf([mfilename ' : ' sprintf('STLN Number of iterations required : %i \n',ite)])
LineBreakLarge()

% if condition(ite) < condition(1) || ite == 1

% Update z and x_ls
vec_z = yy(1:nEntries_z);

% Get z1 and z2
zf = vec_z(1:nCoeffs_fxy);
zg = vec_z(nCoeffs_fxy+1:end);

% Get Coefficients of polynomials f(x,y) and g(x,y)
vec_zf = [zf ; zeros(nZeros_fxy,1)];
vec_zg = [zg ; zeros(nZeros_gxy,1)];

mat_zf = GetAsMatrix(vec_zf,m,m);
mat_zg = GetAsMatrix(vec_zg,n,n);


fxy_lr = fxy + mat_zf;
gxy_lr = gxy + mat_zg;

% Get coefficients of polynomials xu(x,y) and xv(x,y)
% split the vector x into \hat{x}_{1} and \hat{x}_{2}
% Get number of coefficients in x1

% split the vector x into \hat{x}_{1} and \hat{x}_{2}
% Get number of coefficients in x1
nCoefficients_xv = nchoosek(n-k+2,2);
nCoefficients_xu = nchoosek(m-k+2,2);

vec_xvxu = [x_ls(1:idx_col-1) ; -1 ; x_ls(idx_col:end)];

xv = vec_xvxu(1:nCoefficients_xv);
xu = vec_xvxu(nCoefficients_xv + 1 : nCoefficients_xv + nCoefficients_xu);

% Get x1 as a matrix of coefficients for input into BuildT1() function
try
    nZeros_xv = nchoosek(n-k+1,2);
catch
    nZeros_xv = 0;
end

try
    nZeros_xu = nchoosek(m-k+1,2);
catch
    nZeros_xu = 0;
end
% Get vectors of coefficients of v and u
vec_xv = [ xv ; zeros(nZeros_xv,1)];
vec_xu = [ xu ; zeros(nZeros_xu,1)];

% Get as matrix
mat_xv = GetAsMatrix(vec_xv,n-k,n-k);
mat_xu = GetAsMatrix(vec_xu,m-k,m-k);

vxy_lr = mat_xv;
uxy_lr = -1.*mat_xu;
%fxy = fxy ;
%gxy = gxy ;


% else
%
%
% end

plotgraphs_STLN();

end



