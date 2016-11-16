function DPQ = BuildDPQ_STLN(m,n,k,idx_col)
% Build the matrix P_{t} such that any column of the Sylvester matrix
% S_{t}(f,g) is expressed as P_{t}*[f;g]
%
% % Inputs
%
% m : Degree of polynomial f(x,y)
%
% n : Degree of polynomial g(x,y)
%
% k : Index of Sylvester subresultant matrix S_{k}(f,g)
%
% alpha : Optimal value of \alpha
%
% th1 : Optimal value of \theta_{1}
% 
% th2 : Optimal value of \theta_{2}
%
% idx_col : Index of column c_{k} to be removed from S_{k}, and formed by 
% the matrix vector product P*[f;g]
%
% % Outputs 
%
% P : Matrix P such that P*[f;g] = c_{k}

% Get the number of columns in the matrix C_{n-k}(f)
nCols_Tf = nchoosek(n-k+2,2);

% Build the matrix D^{-1}
D = BuildD(m,n-k);

% Build The middle 
if idx_col <= nCols_Tf  % Column is in first partition of S_{k}
    
    % Build G
    G1 = BuildP1_STLN(m,n,k,idx_col);
    G2 = zeros(nchoosek(m+n-k+2,2), nchoosek(n+2,2));
    
    
else    % Opt col is in second partition
    
    % Get the index of the column with respect to the second partition of
    % the Sylvester matrix.
    idx_col = idx_col - nCols_Tf;
    
    % Build G
    G1 = zeros(nchoosek(m+n-k+2,2), nchoosek(m+2,2));
    G2 = BuildP1_STLN(n,m,k,idx_col);
    
end

% Build the matrix Q = [Q1 0 ; 0 Q2]
Q1 = BuildQ1(m);
Q2 = BuildQ1(n);
Q = blkdiag(Q1,Q2);

% Build the matrix P = DGQ
DPQ = D*[G1 G2]*Q;

end