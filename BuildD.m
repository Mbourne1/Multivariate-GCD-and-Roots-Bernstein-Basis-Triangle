function D = BuildD(m,n_t)
% Build the diagonal matrix D^{-1} for the convolution of two polynomials
% f(x,y) and g(x,y) of degrees m and n.

temp_mat = zeros(m+n_t+1,m+n_t+1);

for i = 0:1:m+n_t
    
    for j = 0:1:m+n_t-i
    
        temp_mat(i+1,j+1) = 1./Trinomial(m+n_t,i,j);
    end
    
end


vect = GetAsVector(temp_mat);

% Only want the nchoosek(m+n+2,2) non-zero terms
nNonZeroTerms = nchoosek(m+n_t+2,2);

vect = vect(1:nNonZeroTerms);

D = diag(vect);

end
