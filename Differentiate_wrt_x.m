function [dfxy] = Differentiate_wrt_x(fxy)

[m1,m2] = GetDegree(fxy);
m = m1;

% Get number of non-zero coefficients of f(x,y)
nCoefficients_fxy = nchoosek(m+2,2);

% Get Vector of coefficients
v_fxy = GetAsVector(fxy);
v_fxy = v_fxy(1:nCoefficients_fxy);

% Build matrix of zeros
temp_mat = zeros(m+1,m+1);

nCoefficients_dfxy = nchoosek(m+1,2);

mymat = zeros(nCoefficients_dfxy,nCoefficients_fxy);

count = 1;
for k = 0:1:m-1
    for i = k:-1:0
    
        mat = temp_mat;
        
        j = k-i;
    
        mat(i+1,j+1) = -m;
        mat(i+2,j+1) = m;
    
        temp_vec = GetAsVector(mat)';
        temp_vec = temp_vec(1:nCoefficients_fxy);
        
        mymat(count,:) = temp_vec;
        count = count + 1;
    end
end

% Get nonzero Coefficients of f(x,y)
v_dfxy = mymat * v_fxy;

if m > 1
    nZeros = nchoosek(m,2);
else
    nZeros = 0;
end

% Append zeros to make dfxy fit into a square matrix of size (m-1) x (m-1)
v_dfxy = ...
    [
    v_dfxy;
    zeros(nZeros,1);
    ];


dfxy = GetAsMatrix(v_dfxy,m-1,m-1);


end