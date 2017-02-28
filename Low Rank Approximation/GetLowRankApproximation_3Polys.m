function [fxy_lr, gxy_lr, hxy_lr, uxy_lr, vxy_lr, wxy_lr,...
    alpha_lr, th1_lr, th2_lr] = ...
    GetLowRankApproximation_3Polys(fxy, gxy, hxy, alpha, th1, th2, m, n, o, k)
% Get the low rank approximation of the Sylvester matrix S_{t}(f,g), and
% return the polynomials f_lr(x,y) and g_lr(x,y) which are the perturbed forms of
% input polynomials f(x,y) and g(x,y).
% lr = low rank.
%
% % Inputs.
%
% [fxy, gxy, hxy] : Coefficients of polynomial f(x,y), g(x,y) and h(x,y)
%
% alpha : \alpha
%
% [th1, th2] : \theta_{1}
%
% [m, n, o] : Total degree of f(x,y), g(x,y) and h(x,y)
%
% % Outputs
%
% [fxy_lr, gxy_lr, hxy_lr] : Coefficients of polynomial f_lr(x,y) which is 
% used in the low rank approximation of S_{t}(f,g).
%
% [uxy_lr, vxy_lr, wxy_lr] :
%
% alpha : \alpha
%
% [th1, th2] : \theta_{1}

global SETTINGS

switch SETTINGS.LOW_RANK_APPROX_METHOD
    
    
    case 'None'
        
        % Get f(\omega_{1},\omega_{2})
        fww = GetWithThetas(fxy, m, th1, th2);
        
        % Get g(\omega_{1}, \omega_{2}) 
        gww = GetWithThetas(gxy, n, th1, th2);
        
        % Get h(\omega_{1}, \omega_{2})
        hww = GetWithThetas(hxy, o, th1, th2);
       
        
        
        % Get polynomials u(\omega_{1},\omega_{2}) and
        % v(\omega_{1},\omega_{2}).
        [uww, vww, www] = GetCofactors_3Polys(fww, gww, hww, m, n, o, k);
        
        fxy_lr = fxy;
        gxy_lr = gxy;
        hxy_lr = hxy;
                
        uxy_lr = GetWithoutThetas(uww, m-k, th1, th2);
        vxy_lr = GetWithoutThetas(vww, n-k, th1, th2);
        wxy_lr = GetWithoutThetas(www, o-k, th1, th2);
        
        alpha_lr = alpha;
        th1_lr = th1;
        th2_lr = th2;
        
        SETTINGS.LOW_RANK_APPROX_REQ_ITE = 0;
        
        
    case 'Standard STLN'
        
        error('Not Yet Complete')
        
        fww = GetWithThetas(fxy, m, th1, th2);
        gww = GetWithThetas(gxy, n, th1, th2);
        a_gww = alpha.*gww;
        
        % Get low rank approximation
        [fww_lr, a_gww_lr, uww_lr, vww_lr] = STLN(fww, a_gww,k);
        
        fxy_lr = GetWithoutThetas(fww_lr, m, th1, th2);
        gxy_lr = GetWithoutThetas(a_gww_lr, n, th1, th2) ./ alpha;
        uxy_lr = GetWithoutThetas(uww_lr, m-k, th1, th2);
        vxy_lr = GetWithoutThetas(vww_lr, n-k, th1, th2);
        
        alpha_lr = alpha;
        th1_lr = th1;
        th2_lr = th2;
        
        
        %[uxy_lr,vxy_lr] = GetCofactors(fxy_lr,gxy_lr,t);
        
        switch SETTINGS.PLOT_GRAPHS
            case 'y'
                % Build the Sylvester matrix of f(x,y) and g(x,y)
                S1 = BuildDTQ(fxy, gxy, k);
                S2 = BuildDTQ(fww, a_gww, k);
                S3 = BuildDTQ(fxy_lr, gxy_lr, k);
                S4 = BuildDTQ(fww_lr, a_gww_lr, k);
                
                
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
        
        error('Not Yet Complete')
        
        [fxy_lr, gxy_lr, uxy_lr, vxy_lr, alpha_lr, th1_lr, th2_lr] = ...
            SNTLN(fxy, gxy, alpha, th1, th2, k);
        
        
        
        switch SETTINGS.PLOT_GRAPHS
            case 'y'
                
                % Get f(\omega_{1},\omega_{2}) and g(\omega_{1},\omega_{2})
                fww = GetWithThetas(fxy,m,th1,th2);
                a_gww = alpha .* GetWithThetas(gxy,n,th1,th2);
                
                %
                fww_lr = GetWithThetas(fxy_lr,m,th1_lr,th2_lr);
                a_gww_lr = alpha_lr .*GetWithThetas(gxy_lr,n,th1_lr,th2_lr);
                
                % Build sylvester subresultant matrices
                S1 = BuildDTQ(fxy,gxy,k);
                S2 = BuildDTQ(fww,a_gww,k);
                S3 = BuildDTQ(fxy_lr,gxy_lr,k);
                S4 = BuildDTQ(fww_lr,a_gww_lr,k);
                
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
            case 'n'
        end
end
end