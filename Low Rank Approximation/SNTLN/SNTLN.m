function [fxy_lr,gxy_lr,uxy_lr,vxy_lr] = SNTLN(fxy_matrix,gxy_matrix,alpha,th1,th2,k)
% SNTLN(fxy,gxy,alpha,th1,th2,k)
%
% Compute the low rank approximation of the Sylvester matrix S_{k}(f,g)
% by method of Structured Total Least Norm STLN.
%
% % Inputs.
%
% fxy : Coefficients of polynomial f(x,y)
%
% gxy : Coefficients of polynomial g(x,y)
%
% alpha : Optimal value of \alpha
%
% th1 : Optimal value of \theta_{1}
%
% th2 : Optimal value of \theta_{2}
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
[m,~] = GetDegree(fxy_matrix);

% Get the degree of g(x,y)
[n,~] = GetDegree(gxy_matrix);

% Get number of coefficients in f(x,y)
nCoeff_fxy = nchoosek(m+2,2);

% Get number of coefficients in g(x,y)
nCoeff_gxy = nchoosek(n+2,2);

% Get the number of coefficients in both f(x,y) and g(x,y)
nCoeff_fg = nCoeff_fxy + nCoeff_gxy;

% Get the number of coefficients in v(x,y)
nCoeff_vxy = nchoosek(n-k+2,2);

% Get the number of coefficients in u(x,y)
nCoeff_uxy = nchoosek(m-k+2,2);

% Get the number of coefficients in the unknown vector x, where A_{k}x = c_{k}
nCoeff_x = nCoeff_uxy + nCoeff_vxy - 1;

% Get the number of columns in C_{n-k}(f), the first partition of the
% Sylvester subresultant matrix S_{k}(f,g)
nCols_Tf = nCoeff_vxy;

% Get the number of columns in C_{m-k}(f), the second partititon of the
% Sylvester subresultant matrix S_{k}(f,g)
nCols_Tg = nCoeff_uxy;

% Get the total number of columns in the Sylvester subresultant matrix
% S_{k}(f,g)
nCols_Sk = nCols_Tf + nCols_Tg;

% Create the identity matrix
I = eye(nCols_Sk, nCols_Sk);
M = I;
M(:,idx_col) = [];
e = I(:,idx_col);

%% Preprocessing

% Get the polynomial f(\omega_{1},\omega_{2})
fww_matrix = GetWithThetas(fxy_matrix,m,th1,th2);

% Get the polynomial g(\omega_{1},\omega_{2})
gww_matrix = GetWithThetas(gxy_matrix,n,th1,th2);

%%

% Build the kth Sylvester Subresultant matrix S_{k}(f,g)
% Build the diagonal matrix D^{-1}
D = BuildD(m,n-k);
DTQ_fg = BuildDTQ(fww_matrix,gww_matrix,k1,k2);

% Build the matrix T_{n-k}(f)
T1_fww = BuildT1(fww_matrix,m,n-k);

% Build the matrix T_{m-k}(g)
T1_gww = BuildT1(gww_matrix,n,m-k);

% Build the matrix Q_{k}
Q = BuildQ(m,n,k);

% Build the kth Sylvester subresultant matrix.
DTQ_fg = D*[T1_fww alpha.*T1_gww] * Q;

% Get partial derivatives of f(\omega_{1},\omega_{2}) with respect to \alph
fww_wrt_alpha = zeros(m+1,m+1);

% Get partial derivatives of g(\omega_{1},\omega_{2}) with respect to \alph
alpha_gww_wrt_alpha = gxy_matrix;


% Get index of optimal column for removal from S_{k}(f,g)
idx_col = GetOptimalColumn(DTQ_fg);

% Remove optimal column from S_{k}(f,g)
Ak_fg = DTQ_fg;
Ak_fg(:,idx_col) = [];

% Get c_{k} optimal column removed from S_{k}(f,g)
ck = DTQ_fg(:,idx_col);

% Get least squares solution x_ls
x_ls = SolveAx_b(Ak_fg,ck);

% Get the residual vector associated with LS solution
res_vec = (ck) - (Ak_fg*x_ls);

% Build P_{t} - The matrix P_{t} such that P_{t}z = ht or P_{t}[f;g] = c_{t}
Pk = BuildP(m,n,k,idx_col);

% Obtain the vector \hat{x} which contains x_ls with a zero inserted into
% the position of the optimal column.
vec_xvxu = [x_ls(1:idx_col-1) ; -1 ; x_ls(idx_col:end)];
xk = [x_ls(1:idx_col-1) ; 0 ; x_ls(idx_col:end)];

% Build matrix Y_{k}
Yk = BuildY_STLN(xk,m,n,k);

zf = zeros(nchoosek(m+2,2),1);
zg = zeros(nchoosek(n+2,2),1);
vec_z = [zf;zg];

nZeros_zf = nchoosek(m+1,2);
vec_zf = [zf ; zeros(nZeros_zf,1)];
mat_zf = GetAsMatrix(vec_zf,m,m);

nZeros_zg = nchoosek(n+1,2);
vec_zg = [zg ; zeros(nZeros_zg,1)];
mat_zg = GetAsMatrix(vec_zg,n,n);

T1_zf = BuildT1(mat_zf,m,n-k);
T1_zg = BuildT1(mat_zg,n,m-k);

Sk_zfzg = D*[T1_zf T1_zg] * Q;
Ak_zfzg = Sk_zfzg;
Ak_zfzg(:,idx_col) = [];

% Build matrices H_{z} and H_{x}
H_z = Yk - Pk;
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
    
    
    % Get small petrubations by LSE
    y = LSE(E,f,C,res_vec);
    
    % Increment interation counter
    ite = ite + 1;
       
    % Increment cummulative peturbations
    yy = yy + y;
    
    % Get the entries corresponding to perturbations of f(x,y) and g(x,y)
    delta_zk = y(1:nEntries_z);
    
    % Remove the zk coefficients from the vector y.
    y(1:nCoeff_fxy + nCoeff_gxy) = [];
    
    % Get the entries of y corresponding to x
    delta_xk = y(nEntries_z + 1:end);
    
    % Remove the entries from the vector y
    y(1:nCoeff_x) = [];
    
    % Get the entry in y corresponding to \theta_{1}
    delta_th1 = y(1:1);
    
    % Remove the entry from vector y
    y(1) = [];
    
    % Get the entry in y corresponding to \theta_{2} 
    delta_th2 = y(1:1);
    
    % Remove the entry from vector y
    y(1) = [];
    
    % % Update the variables
    
    % Update variables z_{k}, where z_{k}, where z_{k} are perturbations in
    % the coefficients of f and g.
    zk = zk + delta_zk;
    
    % Update x_{k}
    xk = xk + delta_xk;
    
    % Update alpha
    alpha(ite) = alpha(ite-1) + delta_alpha;
    
    % Update \theta_{1}
    th1(ite) = th1(ite-1) + delta_th1;
    
    % Update \theta_{2}
    th2(ite) = th2(ite-1) + delta_th2;
    
    % %
    % Obtain polynomials in modified Bernstein basis a_{i}\theta^{i} 
    
    %
    fww_matrix = GetWithThetas(fxy_matrix,th1(ite),th2(ite));
    
    %
    gww_matrix = GetWithThetas(gxy_matrix,th1(ite),th2(ite));
    
    %
    DTQ_fg = BuildDTQ(fww_matrix, alpha(ite).*gww_matrix,k1,k2);
    
    % %
    % Get partial derivatives
    
    % Get the partial derivative of f(\omega_{1},\omega_{2}) with respect
    % to \alpha
    fww_wrt_alpha = zeros(m+1,m+1);
    
    % Get the partial derivative of g(\omega_{1},\omega_{2}) with respect
    % to \alpha
    alpha_gww_wrt_alpha = gww_matrix;
    
    % Calculate the partial derivatives of f(\omega_{1},\omega_{2}) with
    % respect to \theta_{1}
    
    % Get the partial derivative of f(w1,w2) with respect to theta_{1}
    fww_wrt_th1 = Differentiate_wrt_theta1(fww_matrix,th1(ite));
    
    % Get the partial derivative of g(w1,w2) with respect to theta_{1}
    gww_wrt_th1 = Differentiate_wrt_theta1(gww_matrix,th1(ite));
        
    % Get the partial derivative of f with respect to theta_{2}
    fww_wrt_th2 = Differentiate_wrt_theta2(fww_matrix,th2(ite));
    
    % Get the partial derivative of g with respect to theta_{2}
    gww_wrt_th2 = Differentiate_wrt_theta2(gww_matrix,th2(ite));
    
    % Calculate the Partial derivative of T with respect to alpha.
    DTQ_wrt_alpha = BuildDTQ(fww_wrt_alpha, alpha_gww_wrt_alpha,k1,k2);
        
    % Calculate the partial derivative of DTQ with respect to theta_{1}
    DTQ_wrt_th1 = BuildDTQ(fww_wrt_th1,alpha(ite).*gww_wrt_th1,k1,k2);
    
    % Calculate the partial derivative of DTQ with respect to theta_{2}
    DTQ_wrt_th2 = BuildDTQ(fww_wrt_th2,alpha(ite).*gww_wrt_th2,k1,k2);
    
    
    % Get z1 and z2
    delta_zf = delta_zk(1:nCoeff_fxy);
    delta_zg = delta_zk(nCoeff_fxy+1:nCoeff_fxy+nCoeff_gxy);
    
    % Update vectors z, z_{f} and z_{g}
    vec_z = vec_z + delta_zk;
    zf =  zf + delta_zf;
    zg =  zg + delta_zg;
    
    % % Build matrix S_{t}(z1,z2)
    vec_zf = [zf ; zeros(nZeros_zf,1)];
    mat_zf = GetAsMatrix(vec_zf,m,m);
    
    vec_zg = [zg ; zeros(nZeros_zg,1)];
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
    x_ls = x_ls + delta_xk;
    
    
    xk = [x_ls(1:idx_col-1) ; 0 ; x_ls(idx_col:end)];
    
    
    % Build the matrix Y_{k}
    Yk = BuildY_STLN(xk,m,n,k);
    
    % % % Build the matrix C in the LSE Problem.
    
    % Build H_z
    H_z = Yk - Pk;
    
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

fprintf([mfilename ' : ' sprintf('STLN Number of iterations required : %i \n',ite)])


% if condition(ite) < condition(1) || ite == 1

% Update z and x_ls
vec_z = yy(1:nEntries_z);

% Get z1 and z2
zf = vec_z(1:nCoeff_fxy);
zg = vec_z(nCoeff_fxy+1:end);

% Get Coefficients of polynomials f(x,y) and g(x,y)
vec_zf = [zf ; zeros(nZeros_zf,1)];
vec_zg = [zg ; zeros(nZeros_zg,1)];

mat_zf = GetAsMatrix(vec_zf,m,m);
mat_zg = GetAsMatrix(vec_zg,n,n);


fxy_lr = fxy_matrix + mat_zf;
gxy_lr = gxy_matrix + mat_zg;

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



