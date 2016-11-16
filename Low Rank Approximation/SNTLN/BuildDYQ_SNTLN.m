function DYQ = BuildDYQ_SNTLN(m, n, k, x, alpha, th1, th2)
% BuildDYQ_SNTLN(m, n, k, x, alpha, th1, th2)
%
% Build the matrix D^{-1}_{m+n-k} * Y_{k}(x1,x2) * Q.
%
% % Inputs
%
% m : Degree of polynomial f(x,y)
%
% n : Degree of polynomial g(x,y)
%
% k : Degree of polynomial d(x,y) and index of kth Sylvester subresultant
% matrix S_{k}(f,g)
%
% x : Vector x
%
% alpha : Optimal value of \alpha
%
% th1 : Optimal value of \theta_{1}
%
% th2 : Optimal value of \theta_{2}
%
% % Outputs
%
% DYQ : Matrix DYQ

% split the vector x into \hat{x}_{1} and \hat{x}_{2}

% Get number of coefficients in x1
nCoeffs_x1 = nchoosek(n-k+2,2);

% Get the number of coefficients in f(x,y)
nCoeffs_fxy = nchoosek(m+2,2);

% Get the number of coefficients in g(x,y)
nCoeffs_gxy = nchoosek(n+2,2);

% Get x1 as a matrix of coefficients for input into BuildT1() function

nZeros_x1 = nchoosek(n-k+1,2);

nZeros_x2 = nchoosek(m-k+1,2);


x1 = x(1:nCoeffs_x1);
x2 = x(nCoeffs_x1 + 1 : end);

% Get vectors of coefficients of x_{v} x_{u} x_{1} and x_{2}
vec_x1 = [ x1 ; zeros(nZeros_x1,1)];
vec_x2 = [ x2 ; zeros(nZeros_x2,1)];

% Get the vectors as matrices of coefficients.
x1_xy = GetAsMatrix(vec_x1,n-k,n-k);
x2_xy = GetAsMatrix(vec_x2,m-k,m-k);

% Build the matrix T_{m}(\hat{x}_{1})
T1_x1 = BuildT1(x1_xy,n-k,m);
T1_x2 = BuildT1(x2_xy,m-k,n);

% % 


% Get a vector with entries \theta_{1}\theta_{2} corresponding to the 
% coefficients of f(\omega_{1},\omega_{2}).
th_f = GetAsVector(GetWithThetas(ones(m+1,m+1),m,th1,th2));
th_f = th_f(1:nCoeffs_fxy);
th_f = diag(th_f);

% Get a vector with entries \theta_{1}\theta_{2} corresponding to the
% coefficients of g(\omega_{1},\omega_{2}).
th_g = GetAsVector(GetWithThetas(ones(n+1,n+1),n,th1,th2));
th_g = th_g(1:nCoeffs_gxy);
th_g = diag(th_g);

% Build a diagonal matrix of thetas
%th_mat = blkdiag(th_f,th_g);

% Build the matrix D_{m+n-k}^{-1}
D = BuildD(m,n-k);

% Build the matrix Q_{m}
Qm = BuildQ1(m);

% Build the matrix Q_{n}
Qn = BuildQ1(n);

% Build DYQ
DYQ = D * [T1_x1*th_f*Qm alpha.*T1_x2*th_g*Qn] ;

end
