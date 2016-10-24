function [hxy] = Bernstein_Multiply(fxy,gxy,m,n)
% Multiply the bivariate bernstein polynomials f(x,y) and g(x,y) of degrees
% m and n respectively.

% Build the matrix D^{-1}T_{1}Q_{1}
D = BuildD(m,n);
T1 = BuildT1(fxy,m,n);
Q1 = BuildQ1(n);

% Get g as a vector
g = GetAsVector(gxy);

% Get number of zero coefficients = \nchoosek(m+1,2)
nNonZeros_gxy = nchoosek(n+2,2);

% Remove the zeros from g
g = g(1:nNonZeros_gxy);

% Comput the vector h, of coefficients of h(x,y)
h = D*T1*Q1 * g;

% Append zeros to hxy to put into a matrix
nZeros_hxy = nchoosek(m+n+1,2);
h = [h; zeros(nZeros_hxy,1)];

% Get the matrix of coefficients of h(x,y)
hxy = GetAsMatrix(h,m+n,m+n);

end