function Q = BuildQ1(m)
% Given the coefficients of polynomial f(x,y), build the diagonal matrix Q
% of trinomial coefficients.

temp_mat = zeros(m+1,m+1);

for i = 0:1:m
    for j = 0:1:m-i
    
        temp_mat(i+1,j+1) = Trinomial(m,i,j);
    
    end
end

vect = GetAsVector(temp_mat);

% Get number of nonzeros
nNonZeros = nchoosek(m+2,2);

%
vect = vect(1:nNonZeros);

Q = diag(vect);


end