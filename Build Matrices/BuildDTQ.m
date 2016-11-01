function S = BuildDTQ(fxy,gxy,k)
% Build the Sylvester matrix D^{-1}T_{k}(f,g)Q where Q is block diagonal
% matrix of matrices Q_{n-k} and Q_{m-k}.

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
S = D* [T1 T2] * Q;

end