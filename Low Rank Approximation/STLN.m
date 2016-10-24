function [fxy,gxy,uxy,vxy] = STLN(fxy,gxy,t)

% Get the degree of f(x,y)
[m,~] = GetDegree(fxy);

% get number of coefficients in f(x,y)
nCoefficients_fxy = nchoosek(m+2,2);

% Get the degree of g(x,y)
[n,~] = GetDegree(gxy);

% Get number of coefficients in g(x,y)
nCoefficients_gxy = nchoosek(n+2,2);

% Build Sylvester Matrix S_{t}(f,g)
D = BuildD(m,n-t);
T1_fx = BuildT1(fxy,m,n-t);
T1_gx = BuildT1(gxy,n,m-t);
Q = BuildQ(m,n,t);
St = D*[T1_fx T1_gx] * Q;

% Get Optimal Column for removal
opt_col = GetOptimalColumn(St);

fprintf([mfilename ' : ' sprintf('Optimal Column for removal : %i \n',opt_col)]);

% Build the matrix M_{q}, the identity matrix of size
nCols_St = nchoosek(m-t+2,2) + nchoosek(n-t+2,2);
Mq = eye(nCols_St);
Mq(:,opt_col) = [];

% Remove optimal column from St(f,g)
At = St;
At(:,opt_col) = [];
ct = St(:,opt_col);

% Get least squares solution x_ls
x_ls = SolveAx_b(At,ct);

% Get the residual vector associated with LS solution
res_vec = ct - (At*x_ls);

% Build P_{t} - The matrix P_{t} such that P_{t}z = ht or P_{t}[f;g] = c_{t}
Pt = BuildPt(m,n,t,opt_col);

% Build the matrix $C_{t}(\hat{x}_{1},\hat{x}_{2})$ such that $C_{t}(x1,x2)*z
% = C_{t}(z1,z2)x$

% Obtain the vector \hat{x} which contains x_ls with a zero inserted into
% the position of the optimal column.
vec_xvxu = [x_ls(1:opt_col-1) ; -1 ; x_ls(opt_col:end)];
vec_x1x2 = [x_ls(1:opt_col-1) ; 0 ; x_ls(opt_col:end)];

% split the vector x into \hat{x}_{1} and \hat{x}_{2}
% Get number of coefficients in x1
nCoefficients_x1 = nchoosek(n-t+2,2);
nCoefficients_x2 = nchoosek(m-t+2,2);

% Get x1 as a matrix of coefficients for input into BuildT1() function
try
    nZeros_xv = nchoosek(n-t+1,2);
catch
    nZeros_xv = 0;
end

try
    nZeros_xu = nchoosek(m-t+1,2);
catch
    nZeros_xu = 0;
end

xv = vec_xvxu(1:nCoefficients_x1);
xu = vec_xvxu(nCoefficients_x1 + 1 : end);
x1 = vec_x1x2(1:nCoefficients_x1);
x2 = vec_x1x2(nCoefficients_x1 + 1 : end);

% Get vectors
vec_xv = [ xv ; zeros(nZeros_xv,1)];
vec_xu = [ xu ; zeros(nZeros_xu,1)];
vec_x1 = [ x1 ; zeros(nZeros_xv,1)];
vec_x2 = [ x2 ; zeros(nZeros_xu,1)];

% Get as matrix
mat_xv = GetAsMatrix(vec_xv,n-t,n-t);
mat_xu = GetAsMatrix(vec_xu,m-t,m-t);
mat_x1 = GetAsMatrix(vec_x1,n-t,n-t);
mat_x2 = GetAsMatrix(vec_x2,m-t,m-t);

% Build the matrix T_{m}(\hat{x}_{1})
T1_x1 = BuildT1(mat_x1,n-t,m);
T1_x2 = BuildT1(mat_x2,m-t,n);


% Build Q_{m}
Qm = BuildQ1(m);

% Build Q_{n}
Qn = BuildQ1(n);


Ct_x1x2 = D * [T1_x1*Qm T1_x2*Qn];


% Build Hz = Y_{t} - C_{t}(x1,x2)
H_z = Ct_x1x2 - Pt;


% %
% %
% %
% Build Hx

% Build matrix C_{t}(f,g)
Ct_fg = D*[T1_fx T1_gx] * Q;

zf = zeros(nchoosek(m+2,2),1);
zg = zeros(nchoosek(n+2,2),1);
vec_zxy = [zf;zg];

% Build matrix C_{t}(z1,z2)
nZeros_zf = nchoosek(m+1,2);
vec_zf = [zf ; zeros(nZeros_zf,1)];
mat_zf = GetAsMatrix(vec_zf,m,m);

nZeros_zg = nchoosek(n+1,2);
vec_zg = [zg ; zeros(nZeros_zg,1)];
mat_zg = GetAsMatrix(vec_zg,n,n);

T1_zf = BuildT1(mat_zf,m,n-t);
T1_zg = BuildT1(mat_zg,n,m-t);

Ct_zfzg = D*[T1_zf T1_zg] * Q;
At_zfzg = Ct_zfzg * Mq;
ht = Ct_zfzg(:,opt_col);
    
% Build Hx
H_x = (Ct_fg + Ct_zfzg) * Mq;

% %
% %
% %
% %
% %
% %
C = [H_z H_x];


% Build the matrix E for LSE Problem
nEntries_x = nchoosek(m-t+2,2) + nchoosek(n-t+2,2) - 1;
nEntries_z = nchoosek(m+2,2) + nchoosek(n+2,2);

E = eye(nEntries_x + nEntries_z);
%E = blkdiag(eye(nEntries_x), zeros(nEntries_z));

% Define the starting vector for the iterations for the LSE problem.
start_point     =   ...
    [...
    vec_zxy;...
    x_ls;
    ];

% Set yy to be the vector which stores all cummulative perturbations.
yy = start_point;

% Set the initial value of vector f to be zero
f = -(yy - start_point);

% Initialise the iteration counter
ite = 1;

% Set the termination criterion
condition(ite) = norm(res_vec)/norm(ct) ;


global SETTINGS
while condition(ite) > SETTINGS.STLN_MAX_ERROR && ite < SETTINGS.STLN_MAX_ITERATIONS
    
    % Increment interation counter
    ite = ite + 1;
    
    % Get small petrubations by LSE
    y_lse = LSE(E,f,C,res_vec);
    
    % Increment cummulative peturbations
    yy = yy + y_lse;

    % update z and x_ls
    delta_zxy = y_lse(1:nEntries_z);
    delta_xls = y_lse(nEntries_z + 1:end);
    
    % Get z1 and z2
    delta_zf = delta_zxy(1:nCoefficients_fxy);
    delta_zg = delta_zxy(nCoefficients_fxy+1:nCoefficients_fxy+nCoefficients_gxy);
    
    vec_zxy = vec_zxy + delta_zxy;
    zf =  zf + delta_zf;
    zg =  zg + delta_zg;
    
    % % Build matrix C_{t}(z1,z2)
    vec_zf = [zf ; zeros(nZeros_zf,1)];
    mat_zf = GetAsMatrix(vec_zf,m,m);

    vec_zg = [zg ; zeros(nZeros_zg,1)];
    mat_zg = GetAsMatrix(vec_zg,n,n);

    T1_zf = BuildT1(mat_zf,m,n-t);
    T1_zg = BuildT1(mat_zg,n,m-t);

    Ct_zfzg = D*[T1_zf T1_zg] * Q;
    
    % Get the matrix At
    At_zfzg = Ct_zfzg * Mq;
    
    % Get the vector ht
    ht = Ct_zfzg(:,opt_col);
    
    % Get x1 and x2
    
    % Obtain the vector \hat{x} which contains x_ls with a zero inserted into
    % the position of the optimal column.
    x_ls = x_ls + delta_xls;
    
    %x_ls = SolveAx_b(At+At_zfzg,ct+ht);
    
    vec_xvxu = [x_ls(1:opt_col-1) ; -1 ; x_ls(opt_col:end)];
    vec_x1x2 = [x_ls(1:opt_col-1) ; 0 ; x_ls(opt_col:end)];
    
    % split the vector x into \hat{x}_{1} and \hat{x}_{2}
    % Get number of coefficients in x1
    xv = vec_xvxu(1:nCoefficients_x1);
    xu = vec_xvxu(nCoefficients_x1 + 1 : nCoefficients_x1 + nCoefficients_x2);
    
    % Get vectors of coefficients of v and u
    vec_xv = [ xv ; zeros(nZeros_xv,1)];
    vec_xu = [ xu ; zeros(nZeros_xu,1)];
    
    % Get as matrix
    mat_xv = GetAsMatrix(vec_xv,n-t,n-t);
    mat_xu = GetAsMatrix(vec_xu,m-t,m-t);
    
    % Build the matrix Y_{t}(x) such that Y_{t}(x) [f;g] = S_{t}(f,g)x    
    x1 = vec_x1x2(1:nCoefficients_x1);
    x2 = vec_x1x2(nCoefficients_x1 + 1 : nCoefficients_x1 + nCoefficients_x2);
    vec_x1 = [x1 ; zeros(nZeros_xv,1)];
    vec_x2 = [x2 ; zeros(nZeros_xu,1)];
    mat_x1 = GetAsMatrix(vec_x1,n-t,n-t);
    mat_x2 = GetAsMatrix(vec_x2,m-t,m-t);
       
    
    % Build the matrix T_{m}(\hat{x}_{1})
    T1_x1 = BuildT1(mat_x1,n-t,m);
    T1_x2 = BuildT1(mat_x2,m-t,n);
    
    
    Ct_x1x2 = D * [T1_x1*Qm T1_x2*Qn];
    
    % % % Build the matrix C in the LSE Problem.
    
    % Build H_z
    H_z = Ct_x1x2 - Pt;
       
    % Build Hx 
    H_x = (Ct_fg + Ct_zfzg) * Mq;
    
    % Build C
    C = [H_z H_x];
    
    % Get the residual vector
    res_vec = (ct + ht) - ((At + At_zfzg)*x_ls);
    
    
    % Update f - used in LSE Problem.
    f = -(yy - start_point);
    
    % Update the termination criterion
    condition(ite) = norm(res_vec)/norm(ct + ht) ;
    
end

fprintf([mfilename ' : ' sprintf('STLN Number of iterations required : %i \n',ite)])


% if condition(ite) < condition(1) || ite == 1
    
    % Update z and x_ls
    vec_zxy = yy(1:nEntries_z);
    
    % Get z1 and z2
    zf = vec_zxy(1:nCoefficients_fxy);
    zg = vec_zxy(nCoefficients_fxy+1:end);
    
    % Build matrix C_{t}(z1,z2)
    vec_zf = [zf ; zeros(nZeros_zf,1)];
    vec_zg = [zg ; zeros(nZeros_zg,1)];
    
    mat_zf = GetAsMatrix(vec_zf,m,m);
    mat_zg = GetAsMatrix(vec_zg,n,n);

    fxy = fxy + mat_zf;
    gxy = gxy + mat_zg;
    
    
    vxy = -1.*mat_xv;
    uxy = mat_xu;
    %fxy = fxy ;
    %gxy = gxy ;
    
    
% else
%     
% 
% end

fig_name = sprintf('%s : condition',mfilename);
figure('name',fig_name)
hold on
plot(log10(condition))
xlabel('iteration');
ylabel('log_{10} condition');
hold off

end


function Pt = BuildPt(m,n,t,opt_col)
%Build the matrix P_{t} such that any column of the Sylvester matrix
%S_{t}(f,g) is expressed as P_{t}*[f;g]

nCols_first_part = nchoosek(n-t+2,2);

if opt_col <= nCols_first_part
    % Opt col is in first partition
    
    % Build D
    D = BuildD(m,n-t);
    
    G1 = BuildG_first_part(m,n,t,opt_col);
    G2 = zeros(nchoosek(m+n-t+2,2), nchoosek(n+2,2));
    
    Q1 = BuildQ1(m);
    Q2 = BuildQ1(n);
    Q = blkdiag(Q1,Q2);
else
    % Opt col is in second partition
    D = BuildD(m,n-t);
    opt_col = opt_col - nCols_first_part;
    G1 = zeros(nchoosek(m+n-t+2,2), nchoosek(m+2,2));
    G2 = BuildG_first_part(n,m,t,opt_col);
    Q1 = BuildQ1(m);
    Q2 = BuildQ1(n);
    Q = blkdiag(Q1,Q2);
    
end

    Pt = D*[G1 G2]*Q;

end

function G = BuildG_first_part(m,n,t,opt_col)
    
    nCols = nchoosek(n-t+2,2);
    try
        nZeros = nchoosek(n-t+1,2);
    catch
        nZeros = 0;
    end
    vec = (1:1:nCols)';
    vec = [vec ; zeros(nZeros,1)];
    
    mat = GetAsMatrix(vec,n-t,n-t);
    
    [row,col] = find(mat==opt_col);
    
    row_index = row - 1;
    col_index = col - 1;
    
    j1 = row_index;
    j2 = col_index;
    
    trinom = Trinomial(n-t,j1,j2);
    
    % Build a matrix of ones corresponding to f(x,y)
    nCoefficients_fxy = nchoosek(m+2,2);
    nZeros_fxy = nchoosek(m+1,2);
    
    fxy = [ones(nCoefficients_fxy,1); zeros(nZeros_fxy,1)];
    fxy = GetAsMatrix(fxy,m,m);
    
    % Get matrix of the product of f(x,y) and v(x,y)
    prod = zeros(m+n-t+1,m+n-t+1);
    
    prod(j1+1:j1+m+1,j2+1:j2+m+1) = fxy;
    
    vec = GetAsVector(prod);
     
    % remove zeros
    nCoefficients_prod = nchoosek(m+n-t+2,2);
    nZeros_prod = nchoosek(m+n-t+1,2);
    
    vec = vec(1:nCoefficients_prod);
    
    A = diag(vec);
    
    % Remove the zero columns
    A(:, find(sum(abs(A)) == 0)) = [] ;
    G = A .* trinom;
end
