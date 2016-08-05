function [maximum,minimum] = GetMaxMin2(a_i1i2,i1,i2,m,n_k)
% Note this function assumes the inclusion of Q in the coefficient matrix.

% Build a 2 dimensional vector to store all occurences of the coefficient
% a_{i_{1},i_{2}} from each of the nchoosek(n-k+2,2) columns in S_{k}(f,g).
mat_ai1i2 = zeros(n_k+1 , n_k+1);

% For each diagonal of v(x,y)
for q = 0:1:n_k
    
    % for each coefficient v_{j1,j2} of v(x,y)
    for j1 = q:-1:0
        
        j2 = q-j1;
        
        entry = a_i1i2 * Trinomial(m,i1,i2) * Trinomial(n_k,j1,j2) ./ Trinomial(m+n_k,i1+j1,i2+j2);
        
        mat_ai1i2(j1+1,j2+1) = entry;
        
    end
end

vec = mat_ai1i2(mat_ai1i2~=0);
mat_ai1i2 = vec;


% take absolute values of A
mat_ai1i2 = abs(mat_ai1i2);

try
[max_r,max_c] = find(mat_ai1i2==max(mat_ai1i2(:)));
[min_r,min_c] = find(mat_ai1i2==min(mat_ai1i2(:)));

% get the maximum and minimum values. Always use (1) since max or min may
% occur more than once, and we are only interested in one of these values.

maximum = mat_ai1i2(max_r(1),max_c(1));
minimum = mat_ai1i2(min_r(1),min_c(1));
catch
    maximum = 0;
    minimum = 0;


end