function [fxy,gxy,dxy,uxy, vxy, t] = o_gcd_mymethod(fxy,gxy,m,n,limits_t)
% Compute the gcd of two input polynomials f(x,y) and g(x,y)
%
%
% % Inputs.
%
% fxy : Coefficients of polynomial f(x,y) in the Bernstein form
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
% fxy :
%
% gxy : 
%
% dxy :
%
% uxy :
%
% vxy :
% 
% t :



% Compute the degree of the GCD
[t,lambda,mu,alpha,th1,th2] = GetGCDDegree(fxy,gxy,m,n);

% % 
% Preprocess
fxy_n = fxy ./ lambda;
gxy_n = gxy ./ mu;

% Get With Thetas
fww = GetWithThetas(fxy_n,m,th1,th2);
gww = GetWithThetas(gxy_n,n,th1,th2);

% Get low rank approximation of S_{t}(f,g) and use perturbed coefficients.
% Update f(w,w) and g(w,w).
[fww,a_gww] = GetLowRankApproximation(fww,alpha.*gww,m,n,t);
gww = a_gww./alpha;


% %
% Get the cofactor polynomials u(x,y) and v(x,y)
[uww,vww] = GetCofactors(fww,alpha.*gww,t);



% %
% Get the GCD polynomial d(x,y)
dww = GetGCDCoefficients(fww,alpha.*gww,uww,vww,m,n,t);

% % 
% Remove thetas
uxy = GetWithoutThetas(uww,m-t,th1,th2);
vxy = GetWithoutThetas(vww,n-t,th1,th2);
dxy = GetWithoutThetas(dww,t,th1,th2);

end