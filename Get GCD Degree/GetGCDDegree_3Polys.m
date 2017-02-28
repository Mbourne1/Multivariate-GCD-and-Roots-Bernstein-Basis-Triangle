function [t, GM_fx, GM_gx, GM_hx, alpha, th1, th2] = ...
    GetGCDDegree_3Polys(fxy, gxy, hxy, m, n, o)
% GetGCDDegree_3Polys(fxy, gxy, hxy, m, n, o)
%
% % Inputs.
%
% [fxy, gxy, hxy] : Coefficients of polynomial f(x,y), g(x,y) and h(x,y)
%
% [m, n, o] : Total degree of f(x,y), g(x,y) and h(x,y)
%
% % Outputs
%
% t :
%
% [GM_fx, GM_gx, GM_hx]
%
% alpha
%
% [th1, th2]

global SETTINGS

vMinimumSingularValues = zeros(1,min([m,n,o]));

for k = 1:1:min([m,n,o])
    
    % %
    % Geometric Mean
    
    % Get Geometric mean of f(x,y) in C_{n-k}(f)
    GM_fx = GetGeometricMean(fxy,m,n-k);
    %GM_fx = 1;
    
    % Get Geometric mean of g(x,y) in C_{m-k}(g)
    GM_gx = GetGeometricMean(gxy,n,m-k);
    %GM_gx = 1;
    
    GM_hx = 1; 
    
    % Divide entries of f(x,y) and g(x,y) by geometric mean
    fxy_n = fxy ./ GM_fx;
    gxy_n = gxy ./ GM_gx;
    hxy_n = hxy ./ GM_hx;
    
    % %
    % Preprocess
    %[alpha,th1,th2] = Preprocess(fxy_n,gxy_n,m,n,k);
    
    alpha = 1;
    th1 = 1;
    th2 = 1;
    
    % Get f(x,y) and g(x,y) with thetas to get f(w_{1},w_{2}) and g(w_{1},w_{2})
    fww = GetWithThetas(fxy_n, m, th1, th2);
    
    gww = GetWithThetas(gxy_n, n, th1, th2);
    
    hww = GetWithThetas(hxy_n, o, th1, th2);
    
    % %
    % Build the Sylvester Matrix
    Sk = BuildSylvesterMatrix_3Polys(fww, gww, hww, m, n, o, k);
    
    % Get vector of singular values
    vSingularValues = svd(Sk);
    
    % Get minimum Singular values
    vMinimumSingularValues(k) = min(vSingularValues);
    
    %
    
    
end

lower_lim = 1;
upper_lim = min(m,n);

if min(m,n) == 1
    t = GetGCDDegree_OneSubresultant(vSingularValues);
else
    t = GetGCDDegree_MultipleSubresultants(vMinimumSingularValues,[lower_lim, upper_lim]);
    
end

switch SETTINGS.PLOT_GRAPHS
    case 'y'
        figure_name = sprintf([mfilename ': ' sprintf('Minimum Singular Values of %s',SETTINGS.SYLVESTER_MATRIX_TYPE)]);
        figure('name',figure_name)
        hold on
        plot(log10(vMinimumSingularValues),'-s')
        xlabel('k')
        ylabel('log_{10} Minimum Singular Values')
        hold off
    case 'n'
end

end