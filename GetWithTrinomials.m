

function fxy_tri = GetWithTrinomials(fxy,m)
% Given the coefficients of the polynomial f(x,y) in Bernstein form, get
% the coefficients in the scaled Bernstein form. That is, coefficients of
% f(x,y) with trinomials included.

for i = 0:1:m
    for j = 0:1:m-i
        
        fxy_tri(i+1,j+1) = fxy(i+1,j+1) .* Trinomial(m,i,j);
        
    end
end

end
