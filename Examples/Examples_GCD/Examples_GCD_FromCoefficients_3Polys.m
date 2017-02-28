function [fxy, gxy, hxy, uxy, vxy, wxy, dxy,...
    m, m1, m2,...
    n, n1, n2,...
    o, o1, o2,...
    t, t1, t2] = ...
    Examples_GCD_FromCoefficients_3Polys(ex_num)
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
% hxy : Coefficients of polynomial h(x,y)
%
% dxy : Coefficients of polynomial d(x,y)
%
% m : Total degree of f(x,y)
%
% n : Total degree of g(x,y)
%
% o : Total degree of h(x,y)
%
% t : Total degree of d(x,y)

syms x y;

addpath('../Examples')

[f_root_mult_arr, g_root_mult_arr, h_root_mult_arr, d_root_mult_arr,...
    u_root_mult_arr, v_root_mult_arr, w_root_mult_arr] = Bivariate_GCD_Examples_3Polys(ex_num);

[fxy] = GetCoefficientsFromSymbolicRoots(f_root_mult_arr);
[gxy] = GetCoefficientsFromSymbolicRoots(g_root_mult_arr);
[hxy] = GetCoefficientsFromSymbolicRoots(h_root_mult_arr);

[dxy] = GetCoefficientsFromSymbolicRoots(d_root_mult_arr);

[uxy] = GetCoefficientsFromSymbolicRoots(u_root_mult_arr);
[vxy] = GetCoefficientsFromSymbolicRoots(v_root_mult_arr);
[wxy] = GetCoefficientsFromSymbolicRoots(w_root_mult_arr);

% Get the symbolic polynomials in power form
symbolic_f = GetSymbolicPoly(f_root_mult_arr);
symbolic_g = GetSymbolicPoly(g_root_mult_arr);
symbolic_h = GetSymbolicPoly(h_root_mult_arr);

symbolic_d = GetSymbolicPoly(d_root_mult_arr);

symbolic_u = GetSymbolicPoly(u_root_mult_arr);
symbolic_v = GetSymbolicPoly(v_root_mult_arr);
symbolic_w = GetSymbolicPoly(w_root_mult_arr);


display(symbolic_f)
display(symbolic_g)
display(symbolic_h)

display(symbolic_d)


% Get the total degree of the polynomials f,g,d when in power form.
[m, m1, m2] = GetDegreeStructure(symbolic_f);
[n, n1, n2] = GetDegreeStructure(symbolic_g);
[o, o1, o2] = GetDegreeStructure(symbolic_h);
[t, t1, t2] = GetDegreeStructure(symbolic_d);


fprintf([mfilename ' : ' sprintf('Total Degree of f(x,y) : %i \n',m)]);
fprintf([mfilename ' : ' sprintf('Total Degree of g(x,y) : %i \n',n)]);
fprintf([mfilename ' : ' sprintf('Total Degree of h(x,y) : %i \n',o)]);
fprintf([mfilename ' : ' sprintf('Total Degree of d(x,y) : %i \n',t)]);

end

function [m,m1,m2] = GetDegreeStructure(symbolic_poly)
%
% % Inputs
%
% symbolic_poly :
%
% % Outputs
%
% m : 
%
% m1 : 
%
% m2 :
%

syms x y

% Get total degree
m = double(feval(symengine, 'degree', symbolic_poly));

% Get relative degree with respect to x
m1 = double(feval(symengine, 'degree', symbolic_poly,x));

% Get relative degree with respect to y
m2 = double(feval(symengine, 'degree', symbolic_poly,y));
end

