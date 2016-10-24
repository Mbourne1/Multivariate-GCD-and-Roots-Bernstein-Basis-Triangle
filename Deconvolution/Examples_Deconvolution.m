function [factors,vMult] = Examples_Deconvolution(ex_num)
% Get the matrix of factors and a vector of corresponding multiplicities.
% Inputs

syms x y;

switch ex_num
    case '1'
        
        % (x + y + 0.5)^7 * (x+y+1.2)^12
        
        factors(1,1) = (x + y + 0.5);
        factors(2,1) = (x + y + 1.2);
        
        % Set the multiplicity of each factor
        vMult = [7; 12];
        
        
    case '2'
        factors(1,1) = (x + y + 0.5);
        factors(2,1) = (x + y + 1.2);
        factors(3,1) = (2*x - y -0.7654);
        
        % Set the multiplicity of each factor
        vMult = [3 ; 7; 11];
        
    case '3'
        % (x + y + 0.5)^7(x + y + 1.2)^10(x + y - 0.15)^15
        factors(1,1) = (x + y + 0.5);
        factors(2,1) = (x + y + 1.2);
        factors(3,1) = (x + y - 0.15);
        
        
        % Set the multiplicity of each factor
        vMult = [7; 10; 15];
        
end
