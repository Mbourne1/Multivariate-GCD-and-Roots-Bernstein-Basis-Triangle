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
% arr_hxy : Array of polynomials h{i}(x,y)

% %
% %
% Form the left hand side matrix

% Get number of polynomials in the array arr_fxy of f_{i}(x,y)
nPolys_arr_fxy = size(arr_fxy,1);

% Get number of polynomials in the array arr_hxy of h_{i}(x,y)
nPolys_arr_hxy = nPolys_arr_fxy - 1;

% Get the degrees of polynomials f_{i}(x,y)
vDeg_arr_fxy = zeros(nPolys_arr_fxy,1);
for i = 1:1:nPolys_arr_fxy
   vDeg_arr_fxy(i)  = GetDegree_Bivariate(arr_fxy{i});
end

vDeg_arr_hxy = vDeg_arr_fxy(1:end-1) - vDeg_arr_fxy(2:end);

% preproc
global SETTINGS
switch SETTINGS.PREPROC_DECONVOLUTIONS
    case 'y'
        %[th1, th2] = GetOptimalTheta(arr_fxy);
        
        th1 = 1;
        th2 = 1;
    case 'n'
        
        th1 = 1;
        th2 = 1;
        
    otherwise 
        
        error([mfilename ' : PREPROC_DECONVOLUTIONS must be y or n'])
end

% Get f(\omega_{1},\omega_{2}) from f(x,y)
arr_fww = cell(nPolys_arr_fxy,1);
for i = 1:1:nPolys_arr_fxy
    
    arr_fww{i} = GetWithThetas(arr_fxy{i},vDeg_arr_fxy(i),th1,th2);
    
end

% Build the matrix C(f1,...,fd)
C_fw = BuildC(arr_fww);


% %
% %
% Form the right hand side vector
vRHS = BuildRHS_vec(arr_fww);


% %
% %
% %
% Get vector of coefficients of the polynomials h_{i}(x,y)
v_hww = SolveAx_b(C_fw,vRHS);

% split solution vector into polynomials h_{i}(x,y)
vDeg_hxy = vDeg_fxy(1:end-1) - vDeg_fxy(2:end);

% Get array of polynomials
arr_hww = Get_hxy(v_hww,vDeg_hxy);

% Get without thetas
arr_hxy = cell(nPolys_arr_hxy,1);
for i = 1:1:nPolys_arr_hxy
    arr_hxy{i} = GetWithoutThetas(arr_hww{i},vDeg_arr_hxy(i),th1,th2);
end


end

function C_fxy = BuildC(arr_fxy)
% Build the matrix C(f1,...,fd)
%
% Inputs.
%
% arr_fxy : Array of polynomials f(x,y)
%
% Outputs.
%
% C_fxy :

% Get number of polynomials in f_{i}(x,y)
nPolys_arr_fxy = size(arr_fxy,1);

% Get degree of each polynomial f_{i}(x,y)
vDeg_fxy = zeros(nPolys_arr_fxy,1);
for i = 1 : 1 : nPolys_arr_fxy
    vDeg_fxy(i) = GetDegree_Bivariate(arr_fxy{i});
end

arr_DT1Q1 = cell(nPolys_arr_fxy, 1);

% For each of the polynomials excluding the first f_{1},...,f_{d}
for i = 2:1:nPolys_arr_fxy
    
    % Get the degree of f{i-1}
    n = vDeg_fxy(i-1);
    
    % Get the degree of f{i}
    m = vDeg_fxy(i);
    
    % Temporarily call the ith entry f(x,y)
    fxy = arr_fxy{i};
    
    % Build the matrix T_{n-m}(f(x,y))
    T1 = BuildT1(fxy,m,n-m);
    
    
    D = BuildD_2Polys(m,n-m);
    Q1 = BuildQ1(n-m);
    
    arr_DT1Q1{i-1} = D*T1*Q1;
end



C_fxy = blkdiag(arr_DT1Q1{:});


end

function vRHS = BuildRHS_vec(arr_fxy)


% Get number of polynomials in array f_{i}(x,y)
nPolys_arr_fxy = size(arr_fxy,1);

% Get the degree of the polynomials f_{i}(x,y)
vDeg_arr_fxy = zeros(nPolys_arr_fxy,1);
for i = 1:1:nPolys_arr_fxy
    vDeg_arr_fxy(i) = GetDegree(arr_fxy{i});
end


% Initialise an array
arr_rhs = cell(nPolys_arr_fxy,1);

% For each polynomial f_{1},...,f_{d} (Note exclude f_{0})
for i = 1:1:nPolys_arr_fxy -1
    
    fxy = GetAsVector(arr_fxy{i});
    
    % Get degree of f(x,y)
    m = vDeg_arr_fxy(i);
    
    % Get number of coefficients in f(x,y)
    nCoefficients_fxy = nchoosek(m+2,2);
    
    % Strip zeros from f(x,y)
    v_fxy = fxy(1:nCoefficients_fxy);
    
    arr_rhs{i} = v_fxy;
    
end

vRHS = cell2mat(arr_rhs);

end

function arr_hxy = Get_hxy(v_hxy,vDeg_arr_hxy)

% Get number of polynomials in array v_{x,y}
nPolys_arr_hxy = size(vDeg_arr_hxy,1);

% Initialise a cell array to store h_{i}(x,y).
arr_hxy = cell(nPolys_arr_hxy,1);


for i = 1:1:nPolys_arr_hxy
    
    % the polynomial h_{i} has degree m_{i-1} - m_{i}
    n = vDeg_arr_hxy(i);
    
    % Number of coefficients in h_{i}(x,y)
    nCoefficients_hxy = nchoosek(n+2,2);
    
    % Vector of coefficients of h_{i}(x,y)
    vec = v_hxy(1:nCoefficients_hxy);
    
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
    v_hxy(1:nCoefficients_hxy) = [];
end


end


