function [fxy,gxy,uxy,vxy,dxy,m,n,t] = Examples_GCD_FromCoefficients(ex_num)
%
% % Inputs
%
% ex_num : Example number
%
%
% % Outputs
%
% fxy : Coefficients of polynomial f(x,y)
%
% gxy : Coefficients of polynomial g(x,y)
%
% dxy : Coefficients of polynomial d(x,y)
%
% m : Total degree of f(x,y)
%
% n : Total degree of g(x,y)
%
% t : Total degree of d(x,y)

syms x y;

addpath('../Examples')
[f,g,d,u,v] = Bivariate_GCD_Examples(ex_num);

% Get the coefficients of the polynomials in Triangular Bernstein form.
fxy = GetCoefficientsFromFactors(f);
gxy = GetCoefficientsFromFactors(g);
dxy = GetCoefficientsFromFactors(d);
uxy = GetCoefficientsFromFactors(u);
vxy = GetCoefficientsFromFactors(v);

% Get the symbolic polynomials in power form
symbolic_f = GetSymbolicPoly(f);
symbolic_g = GetSymbolicPoly(g);
symbolic_d = GetSymbolicPoly(d);

% Get the total degree of the polynomials f,g,d when in power form.
m = double(feval(symengine, 'degree', symbolic_f));
n = double(feval(symengine, 'degree', symbolic_g));
t = double(feval(symengine, 'degree', symbolic_d));

fprintf([mfilename ' : ' sprintf('Total Degree of f(x,y) : %i \n',m)]);
fprintf([mfilename ' : ' sprintf('Total Degree of g(x,y) : %i \n',n)]);
fprintf([mfilename ' : ' sprintf('Total Degree of d(x,y) : %i \n',t)]);

end



