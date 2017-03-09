function [fxy_lr, gxy_lr, uxy_lr, vxy_lr, alpha_lr, th1_lr, th2_lr] = ...
    GetLowRankApproximation_2Polys(fxy, gxy, alpha, th1, th2, m, n, k)
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
% k : Total degree of d(x,y).
%
% % Outputs
%
% fxy_lr : Coefficients of polynomial f_lr(x,y) which is used in the low
% rank approximation of S_{t}(f,g).
%
% gxy_lr : Coefficients of polynomial g_lr(x,y) which is used in the low
% rank approximation of S_{t}(f,g).
%
% uxy_lr :
%
% vxy_lr :
%
% alpha : \alpha
%
% th1 : \theta_{1}
%
% th2 : \theta_{2}

global SETTINGS

switch SETTINGS.LOW_RANK_APPROX_METHOD
    
    
    case 'None'
        
        
        fww = GetWithThetas(fxy, m, th1, th2);
        gww = GetWithThetas(gxy, n, th1, th2);
        a_gww = alpha.*gww;
        
        % Get polynomials u(\omega_{1},\omega_{2}) and
        % v(\omega_{1},\omega_{2}).
        [uww,vww] = GetCofactors_2Polys(fww, a_gww,k);
        
        fxy_lr = fxy;
        gxy_lr = gxy;
        
        uxy_lr = GetWithoutThetas(uww, m-k, th1, th2);
        vxy_lr = GetWithoutThetas(vww, n-k, th1, th2);
        
        alpha_lr = alpha;
        th1_lr = th1;
        th2_lr = th2;
        
        SETTINGS.LOW_RANK_APPROX_REQ_ITE = 0;
        
        
    case 'Standard STLN'
        
        fww = GetWithThetas(fxy, m, th1, th2);
        gww = GetWithThetas(gxy, n, th1, th2);
        a_gww = alpha.*gww;
        
        % Get low rank approximation
        [fww_lr,a_gww_lr,uww_lr,vww_lr] = STLN(fww, a_gww, k);
        
        
        fxy_lr = GetWithoutThetas(fww_lr, m, th1, th2);
        gxy_lr = GetWithoutThetas(a_gww_lr, n, th1, th2) ./ alpha;
        uxy_lr = GetWithoutThetas(uww_lr, m-k, th1, th2);
        vxy_lr = GetWithoutThetas(vww_lr, n-k, th1, th2);
        
        alpha_lr = alpha;
        th1_lr = th1;
        th2_lr = th2;
        
        
        %[uxy_lr,vxy_lr] = GetCofactors(fxy_lr,gxy_lr,t);
        
        if( SETTINGS.PLOT_GRAPHS)
            
            % Build the Sylvester matrix of f(x,y) and g(x,y)
            S1 = BuildDTQ_2Polys(fxy, gxy, k);
            S2 = BuildDTQ_2Polys(fww, a_gww, k);
            S3 = BuildDTQ_2Polys(fxy_lr, gxy_lr, k);
            S4 = BuildDTQ_2Polys(fww_lr, a_gww_lr, k);
            
            
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
            
            
        end
        
    case 'Standard SNTLN'
        
        [fxy_lr, gxy_lr, uxy_lr, vxy_lr, alpha_lr, th1_lr, th2_lr] = ...
            SNTLN(fxy, gxy, alpha, th1, th2, k);
        
        
        
        if( SETTINGS.PLOT_GRAPHS)
            
                
                % Get f(\omega_{1},\omega_{2}) and g(\omega_{1},\omega_{2})
                fww = GetWithThetas(fxy, m, th1, th2);
                a_gww = alpha .* GetWithThetas(gxy, n, th1, th2);
                
                %
                fww_lr = GetWithThetas(fxy_lr, m, th1_lr, th2_lr);
                a_gww_lr = alpha_lr .*GetWithThetas(gxy_lr, n, th1_lr, th2_lr);
                
                % Build sylvester subresultant matrices
                S1 = BuildDTQ_2Polys(fxy, gxy, k);
                S2 = BuildDTQ_2Polys(fww, a_gww, k);
                S3 = BuildDTQ_2Polys(fxy_lr, gxy_lr, k);
                S4 = BuildDTQ_2Polys(fww_lr, a_gww_lr, k);
                
                % Get singular values
                vSingularValues_1 = svd(S1);
                vSingularValues_2 = svd(S2);
                vSingularValues_3 = svd(S3);
                vSingularValues_4 = svd(S4);
                
                % Plot singular values
                figure()
                plot(log10(vSingularValues_1),'-s','DisplayName','f(x,y) g(x,y)')
                hold on
                plot(log10(vSingularValues_2),'-s','DisplayName','f(w,w) g(w,w)')
                plot(log10(vSingularValues_3),'-s','DisplayName','f(x,y)_{lr} g(x,y)_{lr}')
                plot(log10(vSingularValues_4),'-s','DisplayName','f(w,w)_{lr} g(w,w)_{lr}')
                hold off
            
        end
end
end