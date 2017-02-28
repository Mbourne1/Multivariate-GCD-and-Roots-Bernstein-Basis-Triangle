function [fxy, gxy, hxy, uxy, vxy, wxy, dxy,...
    m, m1, m2,...
    n, n1, n2,...
    o, o1, o2,...
    t, t1, t2]  = Examples_GCD_3Polys(ex_num)
% Examples_GCD(ex_num)
% Get the coefficient matrix of the polynomials f(x,y), g(x,y), the GCD
% d(x,y), and cofactors u(x,y) and v(x,y) as well as their corresponding
% degrees.
%
% This file either produces an example from coefficient matrices or 
% 'from roots' where the coefficient matrices are generated.
%
% % Inputs
% 
% ex_num : Example number
%
% % Outputs
%
% [fxy, gxy, hxy] : Coefficient matrix of polynomials f(x,y) g(x,y) and
% h(x,y)
%
% [uxy, vxy, wxy] : Coefficient matrix of polynomials u(x,y) v(x,y) and
% w(x,y)
%
% dxy : Coefficient matrix of polynomial d(x,y)
%
% [m, m1, m2] : Degree structure of polynomial f(x,y)
%
% [n, n1, n2] : Degree structure of polynomial g(x,y)
%
% [o, o1, o2] : Degree structure of polynomial h(x,y)
%
% [t, t1, t2] : Degree structure of polynomial d(x,y)

EXAMPLE_TYPE = 'From Coefficients';

switch EXAMPLE_TYPE
%     case 'From Roots'
%         
%         % Get example polynomials
%         [fxy, gxy,...
%             uxy,vxy,...
%             dxy,...
%             m,m1,m2,...
%             n,n1,n2,...
%             t_exact,t1_exact,t2_exact] = Examples_GCD_FromRoots_3Polys(ex_num);
        
    case 'From Coefficients'
        [fxy, gxy, hxy, ...
            uxy, vxy, wxy,...
            dxy,...
            m, m1, m2,...
            n, n1, n2,...
            o, o1, o2,...
            t, t1, t2] = Examples_GCD_FromCoefficients_3Polys(ex_num);
        
end
end


