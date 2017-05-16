function Sk = BuildSylvesterMatrix_2Polys(fxy, gxy, m, n, k)
% BuildSylvesterMatrix(fxy,gxy,m,n,k)
%
% Build the Sylvester subresultant matrix S_{k}(f,g)
% 
%
% Inputs
%
% fxy : (Matrix) Coefficients of the polynomial f(x,y)
%
% gxy : (Matrix) Coefficients of the polynomial g(x,y)
%
% m : (Int) Total degree of polynomial f(x,y)
%
% n : (Int) Total degree of polynomial g(x,y)
%
% k : (Int) Index of the Sylvester subresultant matrix S_{k} to be constructed.
%
% % Outputs
%
% Sk : (Matrix) kth Sylvester subresultant matrix S_{k}(f,g)

% Build the diagonal matrix D^{-1}
D = BuildD_2Polys(m, n-k);

% Build the matrix T_{n-k}(f(x,y))
T1_fx = BuildT1(fxy, m, n-k);

% Build the matrix T_{m-k}(g(x,y))
T1_gx = BuildT1(gxy, n, m-k);

Q = BuildQ_2Polys(m,n,k);

global SETTINGS
switch SETTINGS.SYLVESTER_BUILD_METHOD
    case 'T'
        
        Sk = [T1_fx T1_gx] ;
        
    case 'DT'
        
        Sk = D*[T1_fx T1_gx];
        
    case 'DTQ'
        
        Sk = D*[T1_fx T1_gx]*Q;
    
    case 'TQ'
    
        Sk = [T1_fx T1_gx]*Q;
        
end


end
