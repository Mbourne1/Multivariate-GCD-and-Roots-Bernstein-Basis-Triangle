
function [fxy, m] = Examples_Roots_FromCoefficients(ex_num)
%
% % Inputs
%
% ex_num : (String)
%
% % Outputs
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% m : Total degree of f(x,y)

syms x y;
addpath(genpath('../Examples'));

f_root_mult_arr = Roots_Examples_Bivariate(ex_num);


[fxy] = GetCoefficientsFromSymbolicRoots(f_root_mult_arr);

% Get the symbolic polynomials in power form
symbolic_f = GetSymbolicPoly(f_root_mult_arr);


display(symbolic_f)

% Get the total degree of the polynomials f,g,d when in power form.
m = double(feval(symengine, 'degree', symbolic_f));

fprintf([mfilename ' : ' sprintf('Total Degree of f(x,y) : %i \n',m)]);


end
