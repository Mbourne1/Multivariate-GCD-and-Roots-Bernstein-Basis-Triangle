function [fxy,gxy,dxy,uxy, vxy, t] = o_gcd_mymethod(fxy,gxy,m,n,limits_t)
% Compute the gcd of two input polynomials f(x,y) and g(x,y), where f(x,y)
% and g(x,y) are bivariate polynomials in Bernstein form, with total degree
% m and n respectively.
%
%
% % Inputs.
%
% fxy : Coefficients of bivariate polynomial f(x,y) in the Bernstein form.
%       
% gxy : Coefficients of polynomial g(x,y) in the Bernstein form
%
% m : Total degree of f(x,y)
%
% n : Total degree of g(x,y)
%
% limits_t :
%
% % Outputs.
%
% fxy : Coefficients of bivariate polynomial f(x,y). 
%
% gxy : Coefficients of bivariate polynomial g(x,y)
%
% dxy : Coefficients of the GCD of f(x,y) and g(x,y)
%
% uxy : Cofactor polynomial u(x,y) where f(x,y)/d(x,y) = u(x,y)
%
% vxy : Cofactor polynomial v(x,y) where g(x,y)/d(x,y) = v(x,y)
% 
% t : Total degree of the GCD d(x,y)

% Note : fxy and gxy are square matrices where entries in the upper left 
% triangle are the coefficients of f(x,y) and g(x,y) respectively. The 
% lower right triangle entries are all zeros.

% Note the coefficients of f(x,y) and g(x,y) on output may differ from the 
% input coefficients to this function since f(x,y) and g(x,y) have been
% normalized by geometric mean, and have been perturbed in the case where
% low rank approximation of the Sylvester matrix S_{t}(f,g) has been used.

% Compute the degree of the GCD
[t,lambda,mu,alpha,th1,th2] = GetGCDDegree(fxy,gxy,m,n);

% % 
% Normalize fxy and gxy to obtain fxy_n and gxy_n
fxy_n = fxy ./ lambda;
gxy_n = gxy ./ mu;

% Get With Thetas
fww = GetWithThetas(fxy_n,m,th1,th2);
gww = GetWithThetas(gxy_n,n,th1,th2);

% Get low rank approximation of S_{t}(f,g) and use perturbed coefficients.
% Update f(w,w) and g(w,w).

[fww,a_gww,uww,vww] = GetLowRankApproximation(fww,alpha.*gww,m,n,t);
gww = a_gww./alpha;

% %
% Get the GCD polynomial d(x,y)
dww = GetGCDCoefficients(fww,alpha.*gww,uww,vww,m,n,t);

% % 
% Remove thetas
uxy = GetWithoutThetas(uww,m-t,th1,th2);
vxy = GetWithoutThetas(vww,n-t,th1,th2);
dxy = GetWithoutThetas(dww,t,th1,th2);
fxy = GetWithoutThetas(fww,m,th1,th2);
gxy = GetWithoutThetas(gww,n,th1,th2);


end