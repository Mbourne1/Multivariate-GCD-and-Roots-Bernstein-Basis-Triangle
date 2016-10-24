
function [fxy,m] = Examples_Roots_FromCoefficients(ex_num)

syms x y;

switch ex_num
    case '1'
        
        f = [
            (x + y -0.273547)   4
            (x + y + 0.6271)    2
            ];
        f = GetFactors(f);
    case '2'
        f = [
            (x - 0.7897)    5
            (y - 0.2323)    5  
            (x - 0.2456)    3
            ];
        f = GetFactors(f);
    case '3'
        f = [
            (x-1.213534)    3
            (x-0.645312)    2
            (y-0.55431)     1
            ];
        f = GetFactors(f);
    case '4'
        f = [
            (x + y -0.273547)   6
            (x + y + 0.6271)    4
            ];
        f = GetFactors(f);
        
    otherwise
        error([mfilename ' : error : Not a valid example number'])
end

% Get the coefficients of f(x,y) in the bernstein basis
fxy = GetCoefficientsFromFactors(f);

% Get the symbolic polynomial in power basis
symbolic_f = GetSymbolicPoly(f);

display(symbolic_f)

% Get the total degree of f(x,y)
m = double(feval(symengine, 'degree', symbolic_f));


fprintf([mfilename ' : ' sprintf('Total Degree of f(x,y) : %i \n',m)]);

end
