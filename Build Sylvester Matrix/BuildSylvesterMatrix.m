function Sk = BuildSylvesterMatrix(fxy,gxy,m,n,k)
% BuildSylvesterMatrix(fxy,gxy,m,n,k)
%
% Build the Sylvester subresultant matrix S_{k}(f,g)
% 
%
% Inputs
%
% fxy : Coefficients of the polynomial f(x,y)
%
% gxy : Coefficients of the polynomial g(x,y)
%
% m : Total degree of polynomial f(x,y)
%
% n : Total degree of polynomial g(x,y)
%
% k : Index of the Sylvester subresultant matrix S_{k} to be constructed.

% Build the diagonal matrix D^{-1}
D = BuildD(m,n-k);

% Build the matrix T_{n-k}(f)
T1_fx = BuildT1(fxy,m,n-k);

% Build the matrix T_{m-k}(g)
T1_gx = BuildT1(gxy,n,m-k);

Q = BuildQ(m,n,k);

global SETTINGS
switch SETTINGS.SYLVESTER_MATRIX_TYPE
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
