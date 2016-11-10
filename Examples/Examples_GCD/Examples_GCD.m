function [fxy,gxy,uxy,vxy,dxy,m,n,t] = Examples_GCD(ex_num)
%
% % Inputs
%
% ex_num : Example Number as a string
%
% % Outputs.
%
% fxy : Coefficients of polynomial f(x,y) 
%
% gxy : Coefficients of polynomial g(x,y)
%
% uxy : Coefficients of polynomial u(x,y)
%
% vxy : Coefficients of polynomial v(x,y)
%
% dxy : Coefficients of polynomial d(x,y)
%
% m : Degree of polynomial f(x,y)
%
% n : Degree of polynomial g(x,y)
%
% t : Degree of polynomial d(x,y)


EXAMPLE_TYPE = 'From Coefficients';

switch EXAMPLE_TYPE
    case 'From Roots'
        [fxy,gxy,uxy,vxy,dxy,m,n,t] = Examples_GCD_FromRoots(ex_num);
        
    case 'From Coefficients'
        [fxy,gxy,uxy,vxy,dxy,m,n,t] = Examples_GCD_FromCoefficients(ex_num);
        
    otherwise
        error('err');
end


end
