
function [symbolic_fxy] = GetSymbolicPoly(arr_f)
% Given the set of symbolic factors of f(x,y) return the symbolic
% polynomial.

symbolic_fxy = arr_f{1};

for i = 2:1:length(arr_f)
    
    symbolic_fxy = symbolic_fxy * arr_f{i};
end


end