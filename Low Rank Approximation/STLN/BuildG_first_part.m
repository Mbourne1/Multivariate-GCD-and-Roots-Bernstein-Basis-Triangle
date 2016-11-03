function G1 = BuildG_first_part(m,n,t,opt_col)
% Build the matrix G1, which forms a partition of the matrix P, and is used
% in the BuildP() function. This function forms part of the STLN package.

% Get the number of columns in the first partition of the Sylvester matrix
nCols_Tf = nchoosek(n-t+2,2);

% Get the number of zeros
try
    nZeros = nchoosek(n-t+1,2);
catch
    nZeros = 0;
end

vec = (1:1:nCols_Tf)';
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
G1 = A .* trinom;
end