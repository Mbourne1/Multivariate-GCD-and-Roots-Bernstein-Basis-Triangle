function [alpha,th1,th2] = Preprocess(fxy,gxy,m,n,k)
% Obtain optimal values of alpha and theta
%
% % Inputs. 
%
% fxy : Coefficients of polynomial f(x,y)
%
% gxy : Coefficients of polynomial g(x,y)
%
% m : Total degree of f(x,y)
%
% n : Total degree of g(x,y)
%
% k : Total degree of d(x,y)

% Get global variables
global SETTINGS

switch SETTINGS.BOOL_ALPHA_THETA
    case 'y'
        
        % Get maximum of each coefficient a_{i,j} from its entries in the Sylvester
        % Matrix S_{k}(f,g)
        [max_fxy,min_fxy] = GetMaxMin(fxy,m,n-k);
        
        % Get maximum of each coefficient b_{i,j} from its entries in the Sylvester
        % Matrix S_{k}(f,g)
        [max_gxy,min_gxy] = GetMaxMin(gxy,n,m-k);
        
        % Get Optimal alpha and theta
        [alpha, th1,th2] = OptimalAlphaTheta(max_fxy,min_fxy, max_gxy,min_gxy,m,n);
        
    case 'n'
        alpha = 1;
        th1 = 1;
        th2 = 1;
    otherwise
        error('err')
end

end