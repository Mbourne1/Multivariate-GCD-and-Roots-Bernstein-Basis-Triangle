function [fxy_lr,gxy_lr,uxy_lr,vxy_lr] = GetLowRankApproximation(fxy,gxy,m,n,t)
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


global SETTINGS

switch SETTINGS.LOW_RANK_APPROX_METHOD
    
    
    case 'None'
        
        fxy_lr = fxy;
        gxy_lr = gxy;
        [uxy_lr,vxy_lr] = GetCofactors(fxy,gxy,t);
        
        return;
        
        
    case 'Standard STLN'

        % Get low rank approximation
        [fxy_lr,gxy_lr,uxy_lr,vxy_lr] = STLN(fxy,gxy,t);
        
        %[uxy_lr,vxy_lr] = GetCofactors(fxy_lr,gxy_lr,t);
               
        switch SETTINGS.PLOT_GRAPHS
            case 'y'
                % Build the Sylvester matrix of f(x,y) and g(x,y)
                D = BuildD(m,n);
                T1 = BuildT1(fxy,m,n);
                T2 = BuildT1(gxy,n,m);
                Q = BuildQ(m,n,0);
                [vSingularValues_1] = svd(D*[T1 T2]*Q);

                % Build the sylvester matrix of f_lr(x,y) and g_lr(x,y)
                T1 = BuildT1(fxy_lr,m,n);
                T2 = BuildT1(gxy_lr,n,m);
                [vSingularValues_2] = svd(D*[T1 T2]*Q);
                
                % Plot the singular values.
                figure('name','STLN')
                plot(log10(vSingularValues_1),'-s','DisplayName','Without STLN');
                hold on
                plot(log10(vSingularValues_2),'-o','displayname','With STLN');
                legend(gca,'show');
                hold off
                
            case 'n'
            otherwise
                error('err');
        end
        
end
end