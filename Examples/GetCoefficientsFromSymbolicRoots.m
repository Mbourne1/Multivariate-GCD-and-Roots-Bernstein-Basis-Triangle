function [fxy] = GetCoefficientsFromSymbolicRoots(root_mult_arr)

% Get the factors of f(x,y)
arr_sym_factors_fxy = GetFactors(root_mult_arr);

% Get the coefficients of f(x,y)
fxy = GetCoefficientsFromFactors(arr_sym_factors_fxy);

end