
m = 2 ;
n_k = 1;


% Build the matrix of i1
mat_i1 = diag(0:1:m+n_k) * ones(m+n_k+1,m+n_k+1);

mat_i2 = ones(m+n_k+1,m+n_k+1) * diag(0:1:m+n_k);

mat_j1 = diag(0:1:n_k) * ones(n_k+1,n_k+1);

mat_j2 = ones(n_k+1,n_k+1) * diag(0:1:n_k);

% Get number of coefficients in f(x,y)
nCoefficients_fxy = nchoosek(m+2,2);
nCoefficients_vxy = nchoosek(n_k+2,2);
nCoefficients_fv = nchoosek(m+n_k+2,2);

v_i1 = GetAsVector(mat_i1);
v_i1 = v_i1(1:nCoefficients_fv);

v_i2 = GetAsVector(mat_i2);
v_i2 = v_i2(1:nCoefficients_fv);

v_j1 = GetAsVector(mat_j1);
v_j1 = v_j1(1:nCoefficients_vxy);

v_j2 = GetAsVector(mat_j2);
v_j2 = v_j2(1:nCoefficients_vxy);

Cf1 = zeros(nCoefficients_fv,nCoefficients_vxy);
Cf2 = zeros(nCoefficients_fv,nCoefficients_vxy);
Cf3 = zeros(nCoefficients_fv,nCoefficients_vxy);
Cf4 = zeros(nCoefficients_fv,nCoefficients_vxy);

temp_mat = zeros(nCoefficients_fv,nCoefficients_vxy);
temp_arr = cell(nCoefficients_fv,nCoefficients_vxy);

for i = 1:1:nCoefficients_fv
    for j = 1:1:nCoefficients_vxy
        Cf1(i,j) = v_i1(i);
        Cf2(i,j) = v_i2(i);
        Cf3(i,j) = v_j1(j);
        Cf4(i,j) = v_j2(j);
        
        i1 = v_i1(i);
        i2 = v_i2(i);
        j1 = v_j1(j);
        j2 = v_j2(j);
        if ( i1-j1 >= 0) && (i1-j1 <= m) && (i2-j2 >= 0) && (i2-j2<=m)
            temp_mat(i,j) = nchoosek(i1,i1-j1);
            
            coef = sprintf('a_{%i,%i}',i1-j1,i2-j2);
            mystr = sprintf('binom{%i}{%i}',i1,i1-j1);
            mystr2 = sprintf('binom{%i}{%i}',i2,i2-j2);
            mystr3 = sprintf('binom{%i}{%i,%i}',m+n_k-i1-i2,m-i1+j1-i2+j2);
            
            temp_arr2{i,j} = sprintf([coef mystr mystr2 mystr3]);
            
        end
    end
end

display(Cf1)
display(Cf2)
display(Cf3)
display(Cf4)
display(temp_mat)
