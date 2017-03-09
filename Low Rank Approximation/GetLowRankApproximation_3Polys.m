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
        
        % Get f(\omega_{1},\omega_{2}), g(\omega_{1},\omega_{2}) and h(\omega_{1},\omega_{2})
        fww = GetWithThetas(fxy, m, th1, th2);
        gww = GetWithThetas(gxy, n, th1, th2);
        hww = GetWithThetas(hxy, o, th1, th2);
               
        % Get polynomials u(\omega_{1},\omega_{2}) and
        % v(\omega_{1},\omega_{2}).
        [uww, vww, www] = GetCofactors_3Polys(fww, gww, hww, m, n, o, k);
        
        fxy_lr = fxy;
        gxy_lr = gxy;
        hxy_lr = hxy;
                
        % Get u(x,y), v(x,y) and w(x,y)
        uxy_lr = GetWithoutThetas(uww, m-k, th1, th2);
        vxy_lr = GetWithoutThetas(vww, n-k, th1, th2);
        wxy_lr = GetWithoutThetas(www, o-k, th1, th2);
        
        alpha_lr = alpha;
        th1_lr = th1;
        th2_lr = th2;
        
        SETTINGS.LOW_RANK_APPROX_REQ_ITE = 0;
        
        
    case 'Standard STLN'
        
        error('Not Yet Complete')
        
       
        
    case 'Standard SNTLN'
        
        error('Not Yet Complete')
       
end
end