

function [] = test()
%
%


m = 5;
n_k = 5;




results_matrix = zeros(m+1, n_k+1);

for i = 0:1: m
    for j = 0:1:n_k
        
        results_matrix(i+1, j+1) = test2(i,j);
        
        trinomial_test(i+1,j+1) = Trinomial(i+j + 1,i,j);
        
    end
end

display(results_matrix);
display(trinomial_test);



fxy = [1 2 3 4 5 6 ; 1 1 1 1 1 0; 1 1 1 1 0 0 ; 1 1 1 0 0 0 ; 1 1 0 0 0 0; 1 0 0 0 0 0 ];
GetArithmeticMean(fxy,m,n_k);

end




function [my_sum] = test2(m,n_k)
%
% % Inputs
%
% m :
%
% n_k :
%
% % Outputs
%
% my_sum

coefficient_sum = zeros(m+1,m+1);

for k = 0:1:m
    for i1 = k:-1:0
        i2 = k - i1;
        
        temp_sum = 0;
        
        for k2 = 0:1:n_k
            
            for j1 = k2:-1:0
                j2 = k2 - j1;
                
                
                temp_sum = temp_sum + ...
                    (...
                    nchoosek(i1 + j1,i1) ...
                    * nchoosek(i2 + j2,i2) ...
                    * nchoosek(m + n_k - i1- i2 - j1 - j2, m - i1 - i2) ...
                    / nchoosek(m + n_k, n_k) ...
                  );
                
                
            end
        end
        
        coefficient_sum(i1+1,i2+1) = temp_sum;
        
    end
end

display(coefficient_sum);
my_sum = temp_sum;

test_sum = (m+n_k+1)*(m+n_k+2) / ((m+1)*(m+2))

end


function lambda = GetArithmeticMean(fxy, m, n_k)
%
% % Inputs
%
% fxy : (Matrix) Coefficients of the polynomial f(x,y)
%
% m : (Int) Degree of f(x,y)
%
% n_k : (Int)
%
% % Outputs
%
% lambda

Tf = BuildDT1Q1(fxy, m, n_k);



lambda1 = ... 
        ( sum(sum(fxy))...
        * (nchoosek(m+n_k+2,2)) ./ ( nchoosek(m+2,2)^2 * nchoosek(n_k+2,2)  ));

vec = Tf(Tf~=0)
lambda2 = mean(Tf(Tf~=0));

lambda1./lambda2




end


function value1 = test3(fxy,m,n_k)




end
