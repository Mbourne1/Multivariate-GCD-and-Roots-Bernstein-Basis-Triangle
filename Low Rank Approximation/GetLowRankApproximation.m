function [fxy_lr,gxy_lr,alpha_lr,th1_lr,th2_lr] = GetLowRankApproximation(fxy,gxy,alpha,th1,th2,m,n,t)
% Get the low rank approximation of the Sylvester matrix S_{t}(f,g), and
% return the polynomials f_lr(x,y) and g_lr(x,y) which are the perturbed forms of
% input polynomials f(x,y) and g(x,y). 
% lr = low rank.
%
% % Inputs.
%
% fxy : Coefficients of polynomial f(x,y).
%
% gxy : Coefficients of polynomial g(x,y).
%
% alpha : \alpha
%
% th1 : \theta_{1}
% 
% th2 : \theta_{2}
%
% m : Total degree of f(x,y).
%
% n : Total degree of g(x,y).
%
% t : Total degree of d(x,y).
%
% % Outputs 
%
% fxy_lr : Coefficients of polynomial f_lr(x,y) which is used in the low
% rank approximation of S_{t}(f,g).
%
% gxy_lr : Coefficients of polynomial g_lr(x,y) which is used in the low
% rank approximation of S_{t}(f,g).
%
% alpha : \alpha
%
% th1 : \theta_{1}
%
% th2 : \theta_{2}

global SETTINGS

switch SETTINGS.LOW_RANK_APPROX_METHOD
    
    
    case 'None'
        
        
        fxy_lr = fxy;
        gxy_lr = gxy;
        alpha_lr = alpha;
        th1_lr = th1;
        th2_lr = th2;
        
        return;
        
        
    case 'Standard STLN'

        fww = GetWithThetas(fxy,m,th1,th2);
        gww = GetWithThetas(gxy,n,th1,th2);
        a_gww = alpha.*gww;
        
        % Get low rank approximation
        [fww_lr,a_gww_lr,uww_lr,vww_lr] = STLN(fww,a_gww,t);
 
        fxy_lr = GetWithoutThetas(fww_lr,m,th1,th2);
        gxy_lr = GetWithoutThetas(a_gww_lr,n,th1,th2) ./ alpha;
        alpha_lr = alpha;
        th1_lr = th1;
        th2_lr = th2;
        
        
        %[uxy_lr,vxy_lr] = GetCofactors(fxy_lr,gxy_lr,t);
               
        switch SETTINGS.PLOT_GRAPHS
            case 'y'
                % Build the Sylvester matrix of f(x,y) and g(x,y)
                S1 = BuildDTQ(fxy,gxy,t);
                S2 = BuildDTQ(fww,a_gww,t);
                S3 = BuildDTQ(fxy_lr,gxy_lr,t);
                S4 = BuildDTQ(fww_lr,a_gww_lr,t);
                
                
                [vSingularValues_1] = svd(S1);
                [vSingularValues_2] = svd(S2);
                [vSingularValues_3] = svd(S3);
                [vSingularValues_4] = svd(S4);
                
                
                % Plot the singular values.
                figure('name','STLN')
                plot(log10(vSingularValues_1),'-s','DisplayName','(fxy,gxy)');
                hold on
                plot(log10(vSingularValues_2),'-s','displayname','(fww,gww)');
                plot(log10(vSingularValues_3),'-s','displayname','(fxy_lr,gxy_lr)');
                plot(log10(vSingularValues_4),'-s','displayname','(fww_lr,gww_lr)');
                legend(gca,'show');
                hold off
                
            case 'n'
            otherwise
                error('err');
        end
        
    case 'Standard SNTLN'
        
        
        
end
end