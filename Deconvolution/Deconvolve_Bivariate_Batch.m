function [arr_hxy] = Deconvolve_Bivariate_Batch(arr_fxy,vDeg_fxy)
% Get the set of polynomials h_{i} given by the deconvolution of the
% polynomials f_{i}, where h_{i} = f_{i-1}/f_{i}
%
% % Inputs.
%
%
% arr_fxy : Array of polynomials f_{i}(x,y)
%
% vDeg_fxy : Vector containing the total degree of the polynomials f{i}
%
%
% % Outputs.
%
%
% arr_hxy :

% %
% %
% Form the left hand side matrix

% Get number of polynomials in the array arr_fxy
nPolys_fxy = size(arr_fxy,1);

nPolys_hxy = nPolys_fxy - 1;

arr_DT1Q1 = cell(nPolys_fxy,1);

% For each of the polynomials excluding the first f_{1},...,f_{d}
for i = 2:1:nPolys_fxy
    
    % Get the degree of f{i-1}
    n = vDeg_fxy(i-1);
    
    % Get the degree of f{i}
    m = vDeg_fxy(i);
    
    % Temporarily call the ith entry f(x,y)
    fxy = arr_fxy{i};
    
    % Build the matrix T_{n-m}(f(x,y))
    T1 = BuildT1(fxy,m,n-m);
    
    
    D = BuildD(m,n-m);
    Q1 = BuildQ1(n-m);
    
    arr_DT1Q1{i-1} = D*T1*Q1;
end



C = blkdiag(arr_DT1Q1{:});


% %
% %
% Form the right hand side vector
arr_rhs = cell(nPolys_fxy,1);

for i = 1:1:nPolys_fxy -1
    
    fxy = GetAsVector(arr_fxy{i});
    
    % Get degree of f(x,y)
    m = vDeg_fxy(i);
    
    % Get number of coefficients in f(x,y)
    nCoefficients_fxy = nchoosek(m+2,2);
    
    % Strip zeros from f(x,y)
    v_fxy = fxy(1:nCoefficients_fxy);
    
    arr_rhs{i} = v_fxy;
    
end

vRHS = cell2mat(arr_rhs);

% %
% %
% %
x_ls = SolveAx_b(C,vRHS);

% split solution vector into polynomials h_{i}(x,y)
vDeg_hxy = vDeg_fxy(1:end-1) - vDeg_fxy(2:end);

arr_hxy = cell(nPolys_hxy,1);
for i = 1:1:nPolys_hxy
    
    % the polynomial h_{i} has degree m_{i-1} - m_{i}
    n = vDeg_hxy(i);
    
    % Number of coefficients in h_{i}(x,y)
    nCoefficients_hxy = nchoosek(n+2,2);
    
    % Vector of coefficients of h_{i}(x,y)
    vec = x_ls(1:nCoefficients_hxy);
    
    % append zeros to vector so we can form a matrix
    try
        nZeros_hxy = nchoosek(n+1,2);
    catch
        nZeros_hxy = 0;
    end
    
    
    vec = ...
        [
        vec;
        zeros(nZeros_hxy,1)
        ];
    
    % Form matrix of coefficients of h_{i}(x,y)
    arr_hxy{i} = GetAsMatrix(vec,n,n);
    
    % Remove coefficients from solution vector x_ls
    x_ls(1:nCoefficients_hxy) = [];
end


end