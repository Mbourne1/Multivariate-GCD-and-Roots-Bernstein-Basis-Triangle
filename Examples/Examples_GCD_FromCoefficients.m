function [fxy,gxy,dxy,m,n,t] = Examples_GCD_FromCoefficients(ex_num)
%
% % Inputs
%
% ex_num : Example number
%
%
% % Outputs
% 
% fxy : Coefficients of polynomial f(x,y)
%
% gxy : Coefficients of polynomial g(x,y)
%
% dxy : Coefficients of polynomial d(x,y)
%
% m : Total degree of f(x,y)
%
% n : Total degree of g(x,y)
%
% t : Total degree of d(x,y)

syms x y;

switch ex_num
    case '1'
        
        u = (y-0.2)^2;
        v = (y-0.4) * (x - 0.3)^3;
        d = (x + 1) * (x + 0.8)^4;
        
        f = u*d;
        g = v*d;
        
    case '2'
        
        u = (x + y + 3.0124);
        v = (x + y + 5.4512);
        d = (x + y + 1) * (x + y + 2) * (x + y + 2);
        
        f = u*d;
        g = v*d;
        
    case '3'
        
        u = (x + y + 3.0124)^(2);
        v = (x + y + 5.4512);
        d = (x + y + 1) * (x + y + 2)^(2) * (x+1)^(3);
        
        f = u*d;
        g = v*d;
        
        
    case '4'
        
        u = (x + y + 0.0124)^(6);
        v = (x + y + 0.4512)^(3);
        d = (x^2 + y^2 + 0.51)^(2) * (x + y + 1.12)^(3) * (x+0.56);
        
        f = u*d;
        g = v*d;
        
    case '5'
        
        d = (x + y + 0.0124)^(6) * (x + y + 0.923547)^2 * (x + y - 0.123456)^4;
        u = (x-1) * (x-2)^2;
        v = (x^2 - 5*x);
        
        f = u*d;
        g = v*d;
        
        
    case '6'
        
        d = (x-0.5) * (x-0.2)^2 * (x-0.3)^3 * (y-0.5)^(6);
        u = (x-0.4445)^4 * (y-0.4)^4;
        v = (x-0.1)^5 * (y-0.2234)^5;
        
        f = u*d;
        g = v*d;
        
    case '7'
        d = (x-0.5) * (x-0.2)^2 * (x-0.3)^3 ;
        u = (x-0.4445)^4 ;
        v = (x-0.1)^5 ;
        
        f = u*d;
        g = v*d;
    otherwise
        error([mfilename ' : error : Not a valid example number'])
end



m = double(feval(symengine, 'degree', f));
n = double(feval(symengine, 'degree', g));
t = double(feval(symengine, 'degree', d));

fxy = double(rot90(coeffs(f,[x,y],'All'),2));
gxy = double(rot90(coeffs(g,[x,y],'All'),2));

uxy = double(rot90(coeffs(u,[x,y],'All'),2));
vxy = double(rot90(coeffs(v,[x,y],'All'),2));
dxy = double(rot90(coeffs(d,[x,y],'All'),2));

[m1,m2] = GetDegree(fxy);
[n1,n2] = GetDegree(gxy);

[m1_t1,m2_t2] = GetDegree(uxy);
[n1_t1,n2_t2] = GetDegree(vxy);
[t1,t2] = GetDegree(dxy);

temp_mat = zeros(m+1,m+1);
temp_mat(1:m1+1,1:m2+1) = fxy;
fxy = temp_mat;

temp_mat = zeros(n+1,n+1);
temp_mat(1:n1+1,1:n2+1) = gxy;
gxy = temp_mat;

temp_mat = zeros(m-t+1,m-t+1);
temp_mat(1:m1_t1+1,1:m2_t2+1) = uxy;
uxy = temp_mat;

temp_mat = zeros(n-t+1,n-t+1);
temp_mat(1:n1_t1+1,1:n2_t2+1) = vxy;
vxy = temp_mat;

temp_mat = zeros(t+1,t+1);
temp_mat(1:t1+1,1:t2+1) = dxy;
dxy = temp_mat;

fxy = PowerToBernstein(fxy,m);
gxy = PowerToBernstein(gxy,n);
uxy = PowerToBernstein(uxy,m-t);
vxy = PowerToBernstein(vxy,n-t);
dxy = PowerToBernstein(dxy,t);

m1 = double(feval(symengine, 'degree', f,x));
m2 = double(feval(symengine, 'degree', f,y));
n1 = double(feval(symengine, 'degree', g,x));
n2 = double(feval(symengine, 'degree', g,y));
t1_exact = double(feval(symengine, 'degree', d,x));
t2_exact = double(feval(symengine, 'degree', d,y));


display(f)
display(g)
display(d)
display(u)
display(v)

fprintf([mfilename ' : ' sprintf('Total Degree of f(x,y) : %i \n',m)]);
fprintf([mfilename ' : ' sprintf('Total Degree of g(x,y) : %i \n',n)]);
fprintf([mfilename ' : ' sprintf('Total Degree of d(x,y) : %i \n',t)]);

end
