
function [fxy,m] = Examples_Roots_FromCoefficients(ex_num)

syms x y;

switch ex_num
    case '1'
        
        f = (x + y -0.273547)^(4) * (x + y + 0.6271)^2;
    
    case '2'
        f = (x-0.7897)^5 * (y-0.2323)^5 * (x-0.2456)^3;
        
    case '3'
        f = (x-1.213534)^3 * (x-0.645312)^2 * (y-0.55431);
    case '4'
        f = (x + y -0.273547)^(6) * (x + y + 0.6271)^4;
    otherwise
        error([mfilename ' : error : Not a valid example number'])
end


% Get the total degree of f(x,y)
m = double(feval(symengine, 'degree', f));

fxy = double(rot90(coeffs(f,[x,y],'All'),2));

% Get the degree of f(x,y) with respect to x and y
[m1,m2] = GetDegree(fxy);

temp_mat = zeros(m+1,m+1);
temp_mat(1:m1+1,1:m2+1) = fxy;
fxy = temp_mat;

fxy = PowerToBernstein(fxy,m);

display(f)


fprintf([mfilename ' : ' sprintf('Total Degree of f(x,y) : %i \n',m)]);

end
