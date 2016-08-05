function D = BuildD(m,n)
% Build the diagonal matrix D^{-1} for the convolution of two polynomials
% f(x,y) and g(x,y) of degrees m and n.

temp_mat = zeros(m+n+1,m+n+1);

for i = 0:1:m+n
    
    for j = 0:1:m+n-i
    
        temp_mat(i+1,j+1) = 1./Trinomial(m+n,i,j);
    end
    
end


vect = GetAsVector(temp_mat);

% Only want the nchoosek(m+n+2,2) non-zero terms
nNonZeroTerms = nchoosek(m+n+2,2);

vect = vect(1:nNonZeroTerms);

D = diag(vect);

end
