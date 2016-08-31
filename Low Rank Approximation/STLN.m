function [fxy,gxy] = STLN(fxy,gxy,t)

% Get the degree of f(x,y)
[m,~] = GetDegree(fxy);

% get number of coefficients in f(x,y)
nCoefficients_fxy = nchoosek(m+2,2);

% Get the degree of g(x,y)
[n,~] = GetDegree(gxy);

% Get number of coefficients in g(x,y)
nCoefficients_gxy = nchoosek(n+2,2);

% Build Sylvester Matrix
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

% Remove optimal column from St
At = St;
At(:,opt_col) = [];
ct = St(:,opt_col);

% Get x_ls
v_x_ls = SolveAx_b(At,ct);

% Set iteration number = 1;
ite = 1;


res_vec = ct - (At*v_x_ls);


% %
% %
% %
% Build Hz

% Build Pt - The matrix P_{t} such that P_{t}z = ht or P_{t}[f;g] = c_{t}
Pt = BuildPt(m,n,t,opt_col);

% Check the BuildYt function
%------------------------------
%ct
test_f = GetAsVector(fxy);
test_f = test_f(1:nCoefficients_fxy);

test_g = GetAsVector(gxy);
test_g = test_g(1:nCoefficients_gxy);

%test_ct = Yt*[test_f;test_g];
%ct./test_ct
% --------------------------

% Build the matrix$C_{t}(\hat{x}_{1},\hat{x}_{2})$ such that $C_{t}(x1,x2)*z
% = C_{t}(z1,z2)x$

% Obtain the vector \hat{x} which contains x_ls with a zero inserted into
% the position of the optimal column.
x = [v_x_ls(1:opt_col-1) ; 0 ; v_x_ls(opt_col:end)];

% split the vector x into \hat{x}_{1} and \hat{x}_{2}
% Get number of coefficients in x1
nCoefficients_x1 = nchoosek(n-t+2,2);

x1 = x(1:nCoefficients_x1);
x2 = x(nCoefficients_x1 + 1 : end);

% Get x1 as a matrix of coefficients for input into BuildT1() function
try
    nZeros_x1 = nchoosek(n-t+1,2);
catch
    nZeros_x1 = 0;
end

try
nZeros_x2 = nchoosek(m-t+1,2);
catch
    nZeros_x2 = 0;
end
x1 = [ x1 ; zeros(nZeros_x1,1)];
x2 = [ x2 ; zeros(nZeros_x2,1)];

% Get as matrix
x1 = GetAsMatrix(x1,n-t,n-t);
x2 = GetAsMatrix(x2,m-t,m-t);

% Build the matrix T_{m}(\hat{x}_{1})
T1_x1 = BuildT1(x1,n-t,m);

% Build the matrix T_{n}(\hat{x}_{2})
T1_x2 = BuildT1(x2,m-t,n);

% Build Q_{m}
Qm = BuildQ1(m);

% Build Q_{n}
Qn = BuildQ1(n);

Ct_x1x2 = D * [T1_x1*Qm T1_x2*Qn];

test1 = Ct_x1x2 * [test_f;test_g];
test2 = At*v_x_ls;
test1./test2;
% Build Hz = Y_{t} - C_{t}(x1,x2)

H_z = Ct_x1x2 - Pt;




% %
% %
% %
% Build Hx

% Build matrix C_{t}(f,g)
Ct_fg = D*[T1_fx T1_gx] * Q;

z1 = zeros(nchoosek(m+2,2),1);
z2 = zeros(nchoosek(n+2,2),1);
v_zxy = [z1;z2];

% Build matrix C_{t}(z1,z2)
nZeros_z1 = nchoosek(m+1,2);
vec_z1 = [z1 ; zeros(nZeros_z1,1)];
mat_z1 = GetAsMatrix(vec_z1,m,m);

nZeros_z2 = nchoosek(n+1,2);
vec_z2 = [z2 ; zeros(nZeros_z2,1)];
mat_z2 = GetAsMatrix(vec_z2,n,n);

T1_z1 = BuildT1(mat_z1,m,n-t);
T1_z2 = BuildT1(mat_z2,n,m-t);

Ct_z1z2 = D*[T1_z1 T1_z2] * Q;
At_z1z2 = Ct_z1z2 * Mq;
ht = Ct_z1z2(:,opt_col);
    
% Build Hx
H_x = (Ct_fg + Ct_z1z2) * Mq;

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
E = eye(nEntries_x+nEntries_z);

% Define the starting vector for the iterations for the LSE problem.
start_point     =   ...
    [...
    v_zxy;...
    v_x_ls;
    ];

% Set yy to be the vector which stores all cummulative perturbations.
yy = start_point;

% Set the initial value of vector p to be zero
f = -(yy - start_point);


% Initialise the iteration counter
ite = 1;

% Set the termination criterion
condition(ite) = norm(res_vec)./norm(ct+ht);


global SETTINGS
while condition(ite) > SETTINGS.STLN_MAX_ERROR && ite < SETTINGS.STLN_MAX_ITERATIONS
    
    % Increment interation counter
    ite = ite + 1;
    
    % Get small petrubations by LSE
    y_lse = LSE(E,f,C,res_vec);
    
    % Increment cummulative peturbations
    yy = yy + y_lse;

    % update z and x_ls
    v_zxy = yy(1:nEntries_z);
    v_x_ls = yy(nEntries_z+1:end);
    
    % Get z1 and z2
    z1 = v_zxy(1:nCoefficients_fxy);
    z2 = v_zxy(nCoefficients_fxy+1:end);
    
    % Build matrix C_{t}(z1,z2)
    vec_z1 = [z1 ; zeros(nZeros_z1,1)];
    mat_z1 = GetAsMatrix(vec_z1,m,m);

    vec_z2 = [z2 ; zeros(nZeros_z2,1)];
    mat_z2 = GetAsMatrix(vec_z2,n,n);

    T1_z1 = BuildT1(mat_z1,m,n-t);
    T1_z2 = BuildT1(mat_z2,n,m-t);

    Ct_z1z2 = D*[T1_z1 T1_z2] * Q;
    
    At_z1z2 = Ct_z1z2 * Mq;
    ht = Ct_z1z2(:,opt_col);
    
    
    % Get x1 and x2
    
    % Obtain the vector \hat{x} which contains x_ls with a zero inserted into
    % the position of the optimal column.
    x = [v_x_ls(1:opt_col-1) ; 0 ; v_x_ls(opt_col:end)];

    % split the vector x into \hat{x}_{1} and \hat{x}_{2}
    % Get number of coefficients in x1

    x1 = x(1:nCoefficients_x1);
    x2 = x(nCoefficients_x1 + 1 : end);
    
    
    % Build the matrix Y_{t}(x) such that Y_{t}(x) [f;g] = S_{t}(f,g)x
    
    vec_x1 = [ x1 ; zeros(nZeros_x1,1)];
    vec_x2 = [ x2 ; zeros(nZeros_x2,1)];

    % Get as matrix
    mat_x1 = GetAsMatrix(vec_x1,n-t,n-t);
    mat_x2 = GetAsMatrix(vec_x2,m-t,m-t);
    
    % Build the matrix T_{m}(\hat{x}_{1})
    T1_x1 = BuildT1(mat_x1,n-t,m);

    % Build the matrix T_{n}(\hat{x}_{2})
    T1_x2 = BuildT1(mat_x2,m-t,n);

    Ct_x1x2 = D * [T1_x1*Qm T1_x2*Qn];
    
    % % % Build the matrix C in the LSE Problem.
    
    % Build H_z
    H_z = Ct_x1x2 - Pt;
       
    % Build Hx 
    H_x = (Ct_fg + Ct_z1z2) * Mq;
    
    % Build C
    C = [H_z H_x];
    
    % Get the residual vector
    res_vec = (ct + ht) - ((At+At_z1z2)*v_x_ls);
    
    
    % Update f - used in LSE Problem.
    f = -(yy - start_point);
    
    % Update the termination criterion
    condition(ite) = norm(res_vec)./norm(ct+ht) ;
    
end

fprintf([mfilename ' : ' sprintf('STLN Number of iterations required : %i \n',ite)])


if condition(ite) < condition(1) || ite == 1
    
    % Update z and x_ls
    v_zxy = yy(1:nEntries_z);
    v_x_ls = yy(nEntries_z+1:end);

    % Get z1 and z2
    z1 = v_zxy(1:nCoefficients_fxy);
    z2 = v_zxy(nCoefficients_fxy+1:end);
    
    % Build matrix C_{t}(z1,z2)
    vec_z1 = [z1 ; zeros(nZeros_z1,1)];
    mat_z1 = GetAsMatrix(vec_z1,m,m);

    vec_z2 = [z2 ; zeros(nZeros_z2,1)];
    mat_z2 = GetAsMatrix(vec_z2,n,n);

    fxy = fxy + mat_z1;
    gxy = gxy + mat_z2;
    
    
    
    
else
    

end


figure()
plot(log10(condition))

end


function Yt = BuildPt(m,n,t,opt_col)

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

    Yt = D*[G1 G2]*Q;

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
