function Sk = BuildSylvesterMatrix_3Polys(fxy, gxy, hxy, m, n, o, k)
% BuildSylvesterMatrix(fxy,gxy,m,n,k)
%
% Build the Sylvester subresultant matrix S_{k}(f,g)
% 
%
% Inputs
%
% [fxy, gxy, hxy] : Coefficients of the polynomial f(x,y), g(x,y) and
% h(x,y)
%
% [m, n, o] : Total degree of polynomial f(x,y), g(x,y) and h(x,y)
%
% k : Index of the Sylvester subresultant matrix S_{k} to be constructed.

% Build the diagonal matrix D^{-1}
D = BuildD_3Polys(m,n,o,k);

% Build the matrix T_{n-k}(f)
T1_fx = BuildT1(fxy,m,n-k);

% Build the matrix T_{m-k}(g)
T1_gx = BuildT1(gxy,n,m-k);

Q = BuildQ_3Polys(m,n,o,k);

global SETTINGS
switch SETTINGS.SYLVESTER_BUILD_METHOD
    case 'T'
        
        
        Sk = [T1_fx T1_gx] ;
        
    case 'DT'
        
        
        % Build the matrix D
        D = BuildD_3Polys(m,n,o,k);
        %D1 = BuildD_2Polys(m,n-k);
        %D2 = BuildD_2Polys(m,o-k);
        
        % Build the matrix T
        T1_f = BuildT1(fxy,m,n-k);
        T2_f = BuildT1(fxy,m,o-k);
        T3_g = BuildT1(gxy,n,m-k);
        T4_h = BuildT1(hxy,o,m-k);
        
        diagonal = blkdiag(T1_f, T2_f);
        column = [T3_g; T4_h];
        T = [diagonal column];

        % Build the matrix DT
        Sk = D*T;
        
    case 'DTQ'
        
        
        Sk = BuildDTQ_3Polys(fxy, gxy, hxy, m, n, o, k);
        
        %Sk = D*[T1_fx T1_gx]*Q;
    
    case 'TQ'
        
        Sk = [T1_fx T1_gx]*Q;
    otherwise
        
end


end
