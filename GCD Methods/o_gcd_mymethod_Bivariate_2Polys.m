function [fxy_o, gxy_o, dxy_o, uxy_o, vxy_o, t] = o_gcd_mymethod_Bivariate_2Polys(fxy, gxy, m, n, limits_t)
% Compute the gcd of two input polynomials f(x,y) and g(x,y), where f(x,y)
% and g(x,y) are bivariate polynomials in Bernstein form, with total degree
% m and n respectively.
%
%
% % Inputs.
%
% [fxy, gxy] : Coefficients of bivariate polynomial f(x,y) in the Bernstein form.
%  
% m : Total degree of f(x,y)
%
% n : Total degree of g(x,y)
%
% limits_t :
%
% % Outputs
%
% [fxy, gxy] : Coefficients of bivariate polynomial f(x,y) and g(x,y)
%
% dxy : Coefficients of the GCD of f(x,y) and g(x,y)
%
% [uxy, vxy] : Cofactor polynomial u(x,y) where f(x,y)/d(x,y) = u(x,y)
% and Cofactor polynomial v(x,y) where g(x,y)/d(x,y) = v(x,y)
% 
% t : Total degree of the GCD d(x,y)

% Note : fxy and gxy are square matrices where entries in the upper left 
% triangle are the coefficients of f(x,y) and g(x,y) respectively. The 
% lower right triangle entries are all zeros.

% Compute the degree of the GCD
[t, GM_fx, GM_gx, alpha, th1, th2] = GetGCDDegree_Bivariate_2Polys(fxy, gxy, m, n, limits_t);

% % 
% Normalize fxy and gxy to obtain fxy_n and gxy_n
fxy_n = fxy ./ GM_fx;
gxy_n = gxy ./ GM_gx;

% Get low rank approximation of S_{t}(f,g) and use perturbed coefficients.
% Update f(w,w) and g(w,w).
[fxy_lr, gxy_lr, uxy_lr, vxy_lr, alpha_lr, th1_lr, th2_lr] = ...
    GetLowRankApproximation_2Polys(fxy_n, gxy_n, alpha, th1, th2, m, n, t);

% %
% Get the GCD polynomial d(x,y)
[fxy_lra, gxy_lra, uxy_lra, vxy_lra, dxy_lra, alpha_lra, th1_lra, th2_lra] = ...
    APF_2Polys(fxy_lr, gxy_lr, uxy_lr, vxy_lr, alpha_lr, th1_lr, th2_lr, m, n, t);

% Get outputs
fxy_o = fxy_lra;
gxy_o = gxy_lra;
dxy_o = dxy_lra;
uxy_o = uxy_lra;
vxy_o = vxy_lra;


end