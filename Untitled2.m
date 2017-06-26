% File used in the development of work for the geometric mean computation

m = 8;
n = 8;
k = 5;

% Bk_version1 
prod_a = 1;

for i1 = 0:1:m
    for i2 = 0:1:m - i1
   
        for j1 = 0:1:n-k
            for j2 = 0:1:n-k-j1
           
                prod_a = prod_a * nchoosek(i1 + j1, j1);
                
            end
        end
        
        
    end
end

power = 1 ./ (nchoosek(m+2,2) * nchoosek(n-k+2,2));
Bk_version1 = prod_a .^ power;


% Bk_version2

prod_a = 1;
for i1 = 0:1:m
    for i2 = 0:1:m-i1
        
        for j1 = 0:1:n-k-1
            for j2 = 0 : 1 : n-k - j1
   
                prod_a = prod_a * nchoosek(i1 + j1, j1);
                
            end
            
        end
    end
end

prod_b = 1;


for i1 = 0:1:m
    for i2 = 0:1:m-i1
   
        prod_b = prod_b * nchoosek(i1 + n - k, n - k);
    
    end
end

Bk_version2 = (prod_a .^ power) * (prod_b .^ power);


% Bk_version3

prod_a = 1;
for i1=0:1:m
    for i2 = 0:1:m - i1
        for j1 = 0:1:n-k-1
            for j2 = 0:1:n-k-1-j1
                prod_a = prod_a * nchoosek(i1+j1,j1);
            end
        end
    end
end

prod_b = 1;
for i1 = 0:1:m
    for i2 = 0:1:m-i1
        for j1 = 0:1:n-k-1
            prod_b = prod_b * nchoosek(i1 + j1, j1);
        end
    end
end

prod_c = 1;
for i1 = 0:1:m
    for i2 = 0:1:m-i1
        prod_c = prod_c * nchoosek(i1 + n - k, n - k);
        
    end
end

Bk_version3 = (prod_a.^ power) * (prod_b.^power) * (prod_c.^power);

% Bk_version4

prod_a = 1;
for i1=0:1:m
    for i2 = 0:1:m - i1
        for j1 = 0:1:n-k-1
            for j2 = 0:1:n-k-1-j1
                prod_a = prod_a * nchoosek(i1+j1,j1);
            end
        end
    end
end

prod_b = 1;
for i1 = 0:1:m
    for i2 = 0:1:m-i1
        for j1 = 0:1:n-k
            prod_b = prod_b * nchoosek(i1 + j1, j1);
        end
    end
end


Bk_version4 = (prod_a.^ power) * (prod_b.^power);


display(Bk_version1)
display(Bk_version2)
display(Bk_version3)
display(Bk_version4)