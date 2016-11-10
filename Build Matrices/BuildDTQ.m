function Sk = BuildDTQ(fxy,gxy,k)
% BuildDTQ(fxy,gxy,k)
%
% Build the Sylvester matrix D^{-1}T_{k}(f,g)Q where Q is block diagonal
% matrix of matrices Q_{n-k} and Q_{m-k}.
%
% % Inputs
%
% fxy : Coefficients of polynomial f(x,y)
%
% gxy : Coefficients of polynomial g(x,y)
%
% k : index k, determines which Sylvester subresultant matrix is being
% constructed. S_{k}(f,g)
%
% % Outputs
%
% Sk : Sylvester subresultant matrix S_{k}(f,g)

% Get degree of f(x,y) and g(x,y)
m = GetDegree(fxy);
n = GetDegree(gxy);

% Build the matrix D^{-1}
D = BuildD(m,n-k);

% Build the matrices T_{n-k}(f) and T_{m-k}(g)
T1 = BuildT1(fxy,m,n-k);
T2 = BuildT1(gxy,n,m-k);

% Build the block diagonal matrix Q
Q1 = BuildQ1(n-k);
Q2 = BuildQ1(m-k);
Q = blkdiag(Q1,Q2);

% Construct the k-th Sylvester subresultant matrix S_{k}(f,g)
Sk = D* [T1 T2] * Q;

end