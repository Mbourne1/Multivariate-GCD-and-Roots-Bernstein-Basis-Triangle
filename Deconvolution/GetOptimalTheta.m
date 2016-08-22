function [th1,th2] = GetOptimalTheta(arr_fxy)

f = [1 -1 0 0];

% Get the number of polynomials in the array 
nPolys_arr_fxy = size(arr_fxy,1);

% Get the degree of the array of polynomials
vDeg_arr_fxy = zeros(nPolys_arr_fxy,1);
for i = 1:1:nPolys_arr_fxy
    vDeg_arr_fxy(i) = GetDegree(arr_fxy{i});
end

vDeg_arr_hxy = vDeg_arr_fxy(1:end-1) - vDeg_arr_fxy(2:end);

% Create an array which stores the maximum of each coefficient from each
% polynomial f_{i}(x,y)
arr_max_fxy = cell(nPolys_arr_fxy,1);
arr_min_fxy = cell(nPolys_arr_fxy,1);

for i = 2:1:nPolys_arr_fxy
    [arr_max_fxy{i,1} , arr_min_fxy{i,1} ] = GetMaxMin_fxy(arr_fxy{i},vDeg_arr_hxy(i-1));
end

% Initialise cell arrays rho and tau
rho = cell(nPolys_arr_fxy,1);
tau = cell(nPolys_arr_fxy,1);

for i = 2:1:nPolys_arr_fxy
    
    m = GetDegree(arr_fxy{i});
    nCoefficients_fxy = nchoosek(m+2,2);
    
    temp_max = GetAsVector(arr_max_fxy{i});
    temp_max = temp_max(1:nCoefficients_fxy);
    
    temp_min = GetAsVector(arr_min_fxy{i});
    temp_min = temp_min(1:nCoefficients_fxy);
    
    rho{i} = temp_max;
    tau{i} = temp_min;
end

rho_vec = cell2mat(rho);
tau_vec = cell2mat(tau);

% Create an array which stores the minimum of each coefficient from each
% polynomial f_{i}(x,y)
arr_min_fxy = cell(nPolys_arr_fxy,1);


% for each of the polynomials f_{1}...,f_{d}
arr_A = cell(nPolys_arr_fxy,1);
arr_B = cell(nPolys_arr_fxy,1);
for i = 2:1:nPolys_arr_fxy
    
    % Get degree of f_{i}(x)
    m = vDeg_arr_fxy(i);
    
    % Get number of coefficients in f_{i}(x)
    nCoefficients_fxy = nchoosek(m+2,2);
    
    % Build part of the matrix
    arr_A{i,1} = [ones(nCoefficients_fxy,1) zeros(nCoefficients_fxy,1) -1.*GetPairs(m) ];
    arr_B{i,1} = [zeros(nCoefficients_fxy,1) -1.*ones(nCoefficients_fxy,1) 1.*GetPairs(m)];
end

A = cell2mat(arr_A);
B = cell2mat(arr_B);

LHS_Matrix = [A;B];

RHS_Vec = [rho_vec ; - tau_vec];

x = linprog(f,-LHS_Matrix,-RHS_Vec);


try
    th1 = 10^x(3);
    th2 = 10^x(4);
catch
    fprintf('Failed to optimize\n')
    th1 = 1;
    th2 = 1;
end

end

function v_pairs = GetPairs(m)

nCoefficients_fxy = nchoosek(m+2,2);

% Get index for i1
mat_rows = diag(0:1:m) * fliplr(triu(ones(m+1,m+1),0));
v_rows = GetAsVector(mat_rows);
v_rows = v_rows(1:nCoefficients_fxy);

% Get index for i2
mat_cols = fliplr(triu(ones(m+1,m+1),0)) *diag(0:1:m) ;
v_cols = GetAsVector(mat_cols);
v_cols = v_cols(1:nCoefficients_fxy);

v_pairs = [v_rows v_cols];

end

function [max_matrix,min_matrix] = GetMaxMin_fxy(fxy,n)
% Get maximum occurence of each coefficient 

% Get the degree of f(x,y)
m = GetDegree(fxy);

max_matrix = zeros(m+1,m+1);
min_matrix = zeros(m+1,m+1);
for i = 0:1:m
    for i1 = i:-1:0

        i2 = i - i1;
        
        % Get the coefficient
        a_i1i2 = fxy(i1+1,i2+1);
        num_binom_1 = Trinomial(m,i1,i2);
        
        
        ai1i2 = zeros(n+1,n+1);
        for j = 0:1:n
            for j1 = j:-1:0
                j2 = j - j1;
                num_binom_2 = Trinomial(n,j1,j2);
                den_binom_1 = Trinomial(m+n,i1+j1,i2+j2);
                
                ai1i2(j1+1,j2+1) = log10(abs(a_i1i2 * num_binom_1 * num_binom_2 ./ den_binom_1));
                
            end
        end
        
        
        % Get As Vector and remove zeros
        nCoefficients_hxy = nchoosek(n+2,2);
        v_ai1i2 = GetAsVector(ai1i2);
        v_ai1i2 = v_ai1i2(1:nCoefficients_hxy);
        
        max_ai1i2 = max(v_ai1i2);
        min_ai1i2 = min(v_ai1i2);
        
        max_matrix(i1+1,i2+1) = max_ai1i2;
        min_matrix(i1+1,i2+1) = min_ai1i2;
        
    end
end






end