function [t,lambda,mu,alpha,th1,th2] = GetGCDDegree(fxy,gxy,m,n)
% GetGCDDegree(fxy,gxy,m,n)
%
% Inputs.
%
% fxy : Coefficients of polynomial f(x,y)
%
% gxy : Coefficients of polynomial g(x,y)
%
% m : Total degree of f(x,y)
%
% n : Total degree of g(x,y)

global SETTINGS

vMinimumSingularValues = zeros(1,min(m,n));

for k = 1:1:min(m,n)
    
    % %
    % Geometric Mean
    
    % Get Geometric mean of f(x,y) in C_{n-k}(f)
    lambda = GetGeometricMean(fxy,m,n-k);
    
    % Get Geometric mean of g(x,y) in C_{m-k}(g)
    mu = GetGeometricMean(gxy,n,m-k);
    
    % Divide entries of f(x,y) and g(x,y) by geometric mean
    fxy_n = fxy ./ lambda;
    gxy_n = gxy ./ mu;
    
    % %
    % Preprocess
    [alpha,th1,th2] = Preprocess(fxy_n,gxy_n,m,n,k);
    
    % Get f(x,y) and g(x,y) with thetas to get f(w_{1},w_{2}) and g(w_{1},w_{2})
    fww = GetWithThetas(fxy_n,m,th1,th2);
    
    gww = GetWithThetas(gxy_n,n,th1,th2);
    
    % %
    % Build the Sylvester Matrix
    
    % Build the diagonal matrix D^{-1}
    D = BuildD(m,n-k);
    
    % Build the matrix T_{n-k}(f)
    T1 = BuildT1(fww,m,n-k);
    
    % Build the matrix T_{m-k}(g)
    T2 = BuildT1(alpha.*gww,n,m-k);
    
    % Build the matrix Q
    Q = BuildQ(m,n,k);
    
    
    % Build the Sylvester matrix S_{k}
    Sk = D * [T1 T2] * Q;
    
    % Get vector of singular values
    vSingularValues = svd(Sk);
    
    % Get minimum Singular values
    vMinimumSingularValues(k) = min(vSingularValues);
    
    %
    if(k == 1)
        
        switch SETTINGS.PLOT_GRAPHS
            case 'y'
                figure('name',sprintf([mfilename 'Singular Values S_{1}']));
                plot(log10(svd(Sk)),'-s');
                hold off
            case 'n'
            otherwise
                error('err')
                
        end
        
    end
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
        figure_name = sprintf([mfilename ': ' 'Minimum Singular Values']);
        figure('name',figure_name)
        hold on
        plot(log10(vMinimumSingularValues),'-s')
        xlabel('k')
        ylabel('log_{10} Minimum Singular Values')
        hold off
    case 'n'
end

end