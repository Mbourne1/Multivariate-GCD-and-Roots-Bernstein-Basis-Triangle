function [] = BernsteinMultiply()
ex = 2
switch ex
    case 1
        syms x y;
        
        f = (x + 11)^2 * (y+2) ;
        g = (x + 11)^2 * (y+2) ;
        
        % Get total degree of f(x,y) and g(x,y)
        m = double(feval(symengine, 'degree', f));
        n = double(feval(symengine, 'degree', g));
        
        % Get coefficients of f(x,y) and g(x,y) in power basis.
        fxy = double(rot90(coeffs(f,[x,y],'All'),2));
        gxy = double(rot90(coeffs(g,[x,y],'All'),2));
        
        % Get degree of f(x,y)
        [m1,m2] = GetDegree(fxy);
        
        % Get degree of g(x,y)
        [n1,n2] = GetDegree(gxy);
        
        % Put f(x,y) and g(x,y) into matrices of dimension m+1 x m+1 and n+1 x n+1
        temp_mat = zeros(m+1,m+1);
        temp_mat(1:m1+1,1:m2+1) = fxy;
        fxy = temp_mat;
        
        temp_mat = zeros(n+1,n+1);
        temp_mat(1:n1+1,1:n2+1) = gxy;
        gxy = temp_mat;
        
        % Get coefficients of f(x,y) and g(x,y) in Bernstein basis.
        fxy = PowerToBernstein(fxy,m);
        gxy = PowerToBernstein(gxy,n);
    case 2
        
        fxy = ...
            [
            3 2 2
            2 5 0
            1 0 0
            ];
        gxy = ...
            [
            1 2
            3 0
            ];
        
        m = 2;
        n = 1;
        
end

% Build D^{-1}
D = BuildD(m,n);

% Build C(f)
Cf = BuildC(fxy,m,n);

% Build Q
Q = BuildQ1(n);

% Get the vector of coefficients of g(x,y)
v_gxy = GetAsVector(gxy);
nCoefficients_gxy = nchoosek(n+2,2);
v_gxy = v_gxy(1:nCoefficients_gxy);

display(D)
display(Cf)
display(Q)
display(v_gxy)

D*Cf*Q*v_gxy

end