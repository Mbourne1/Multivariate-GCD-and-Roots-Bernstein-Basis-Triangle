% This experiment looks at the scaling of the coefficients of the
% polynomial f(x) in the Sylvester matrix S_{k}(f,g)

function [] = CoefficientScaling(m,n,k)



nNonZeros_fxy = nchoosek(m + 2, 2);
nZeros_fxy = nchoosek(m + 1, 2);
fxy_vec = [ones(nNonZeros_fxy, 1); zeros(nZeros_fxy,1)];
fxy = GetAsMatrix(fxy_vec, m, m);

arrT = cell(m + 1, m + 1);
arrDT = cell(m + 1, m + 1);
arrDTQ = cell(m + 1, m + 1);
arrTQ = cell(m + 1, m + 1);
% Scaling Effect of T

for k1 = 0 : 1  : m
    
    for i1 = k1 : -1 : 0
        
        i2 = k1 - i1;
        
        
        for k2 = 0:1:n-k
            
            for j1 = k2: -1 : 0
                
                j2 = k2 - j1;
                
                
                arrT{i1+1, i2+1}(j1 + 1, j2 + 1) = Trinomial(m, i1, i2);
                
                arrDT{i1+1, i2+1}(j1 + 1, j2 + 1) = ...
                    Trinomial(m, i1, i2)...
                    ./ ...
                    Trinomial(m + n - k, i1 + j1, i2 + j2);
                
                arrTQ{i1+1, i2+1}(j1 + 1, j2 + 1) = ...
                    Trinomial(m, i1, i2) ...
                    * ...
                    Trinomial(n - k, j1, j2);
                
                arrDTQ{i1+1, i2+1}(j1 + 1, j2 + 1) = ...
                    Trinomial(m, i1, i2) ...
                    * ...
                    Trinomial(n - k, j1, j2) ...
                    ./ ...
                    Trinomial(m + n - k, i1 + j1, i2 + j2);
                
                arrDTQ_DenominatorsRemoved{i1+1, i2+1}(j1 + 1, j2 + 1) = ...
                    nchoosek(i1 + j1, i1) ...
                    * ...
                    nchoosek(i2 + j2, i2) ...
                    / ... 
                    nchoosek(m + n - k - i1 - i2 - j1 - j2, m - i1 - i2) ...
                    .* ...
                    nchoosek(m + n - k, n - k);
               
                
                
                
            end
        end
        
        
        
        
        
    end
    
    
end

PlotValues(arrT, m, n, k);
PlotValues(arrDT, m, n, k);
PlotValues(arrTQ, m , n, k);
PlotValues(arrDTQ, m, n, k);
PlotValues(arrDTQ_DenominatorsRemoved, m, n, k);


end




function [] = PlotValues(arrT, m, n, k)

figure()

nNonZeros_fxy = nchoosek(m + 2, 2);
nNonZeros_vxy = nchoosek(n - k + 2, 2);
nZeros_fxy = nchoosek(m + 1, 2);




hold on
for k1 = 0 : 1  : m
    
    for i1 = k1 : -1 : 0
        
        i2 = k1 - i1;
        
        temp_vec = GetAsVector(arrT{i1+1, i2+1});
        temp_vec = temp_vec(1: nNonZeros_vxy);
        plot((temp_vec));
        
    end
end

hold off


end



