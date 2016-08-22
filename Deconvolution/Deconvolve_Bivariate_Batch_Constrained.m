function [arr_hxy] = Deconvolve_Bivariate_Batch_Constrained(arr_fxy,vDegt_fxy)
%
%
% % Inputs
%
% arr_fxy : Array of polynomials f_{i}(x,y) in Bernstein form.
%
% vDegt_fxy : Vector of total degrees of the polynomials f_{i}(x,y).


vDeg_arr_fxy = vDegt_fxy;
vDeg_arr_hxy = diff(vDeg_arr_fxy);
vDeg_arr_wxy = diff([vDeg_arr_hxy; 0]);
vMult = find(vDeg_arr_wxy~=0);
unique_vMult = unique(vMult);

% Get the number of polynomials in arr_fxy of f_{i}(x,y)
nPolys_arr_fxy = size(arr_fxy,1);

% Get the number of polynomials in the array of h_{i}(x,y)
nPolys_arr_hxy = nPolys_arr_fxy - 1;

vDeg_arr_fxy = zeros(nPolys_arr_fxy,1);
for i = 1:1:nPolys_arr_fxy
   vDeg_arr_fxy(i)  = GetDegree(arr_fxy{i});
end

vDeg_arr_hxy = vDeg_arr_fxy(1:end-1) - vDeg_arr_fxy(2:end);


% preproc
global SETTINGS
switch SETTINGS.PREPROC_DECONVOLUTIONS
    case 'y'
        [th1, th2] = GetOptimalTheta(arr_fxy);
    case 'n'
        th1 = 1;
        th2 = 1;
    otherwise
end

% Get f(\omega_{1},\omega_{2}) from f(x,y)
arr_fww = cell(nPolys_arr_fxy,1);
for i = 1:1:nPolys_arr_fxy
    m = vDeg_arr_fxy(i);
    arr_fww{i} = GetWithThetas(arr_fxy{i},m,th1,th2);
end


% %
% %
% Build the LHS Matrix
DT_fwwQ = BuildDTQ(arr_fww,vDeg_arr_fxy);

% %
% %
% Build the RHS Vector
rhs_fww = BuildRHS_vec(arr_fww);


x_ls = SolveAx_b(DT_fwwQ,rhs_fww);

unique_vMult = unique(vMult);

arr_pww = cell(length(unique_vMult),1);

for i = 1:1:length(unique_vMult)
    
    
    mult = unique_vMult(i);
    
    % Get degree of p(x,y)
    deg_px = vDegt_fxy(mult) - vDegt_fxy(mult+1);
    
    % Get number of coefficients in p(x,y)
    nCoefficients_px = nchoosek(deg_px+2,2);
    
    % Get coefficients of p(x,y) from x_ls
    vec_px = x_ls(1:nCoefficients_px);
    
    % Remove coefficients
    x_ls(1:nCoefficients_px) =[];
    
    
    nZeros_px = nchoosek(deg_px+1,2);
    
    vec_px = ...
        [
            vec_px;
            zeros(nZeros_px,1)
        ];
    
    arr_pww{i,1} = GetAsMatrix(vec_px,deg_px,deg_px);
    
    
    
end



% Get number of polynomials in the array
nPolys_arr_pxy = size(arr_pww,1);

% Get array
arr_hww = Get_hxy(arr_pww,unique_vMult);


% Get degree of polynomials h(x,y)
vDeg_arr_hxy = zeros(nPolys_arr_hxy,1);

for i = 1:1:nPolys_arr_hxy
    vDeg_arr_hxy(i) = GetDegree(arr_hww{i});
end



% % Get without thetas

arr_hxy = cell(nPolys_arr_hxy,1);
for i = 1:1:nPolys_arr_hxy
    
    % Get degree of the polynomial
    m = vDeg_arr_hxy(i);
    
    % Get h_{i}(x,y)
    arr_hxy{i} = GetWithoutThetas(arr_hww{i},m,th1,th2);
end

end


function LHS_Matrix = BuildDTQ(arr_fxy,vDegt_fxy)

vDeg_arr_fxy = vDegt_fxy;
vDeg_arr_hxy = diff(vDeg_arr_fxy);
vDeg_wxy = diff([vDeg_arr_hxy; 0]);
vMult = find(vDeg_wxy~=0);


% Get number of distinct polynomials h_{i}(x)
nDistinct_hx = length(vMult);

for i = 1:1:nDistinct_hx
    
    if i>1
        old_mult = vMult(i-1);
    else
        old_mult = 0;
    end
    
    new_mult = vMult(i);
    
    Cf{i} = [];
    
    for j = (old_mult+1+1) : 1 : (new_mult+1)
        
        % Get coefficients of previous f(x,y)
        fxy_prev = arr_fxy{j-1};
        
        % Get the degree of previous f(x,y)
        deg_fx_prev = vDeg_arr_fxy(j-1);
        
        % Get polynomial f(x,y)
        fxy = arr_fxy{j};
        
        % Get the degree of f(x,y)
        m = vDeg_arr_fxy(j);
        
        % Get the degree of polynomial h_{i}
        deg_hx = deg_fx_prev - m;
        
        % Build the cauchy like matrix
        Tf{j} = BuildT1(fxy,m,deg_hx);
        
        % Build the diagonal matrix D
        D{j} = BuildD(m,deg_hx);
        
        % Stack beneath all other T_f
        Cf{i} = [Cf{i} ; D{j}*Tf{j}];
        
    end
    
    Q{i} = BuildQ1(deg_hx);
    
    DTQ{i} = Cf{i} * Q{i};
end

LHS_Matrix = blkdiag(DTQ{:});

end

function rhs_fxy = BuildRHS_vec(arr_fxy)

% Get number of polynomials in array of f_{i}(x,y)
nPolys_arr_fxy = size(arr_fxy,1);

% Get the degree of the polynomials f_{i}(x,y)
vDeg_arr_fxy = zeros(nPolys_arr_fxy,1);
for i = 1:1:nPolys_arr_fxy 
    vDeg_arr_fxy(i) = GetDegree(arr_fxy{i});
end


% Initialise cell array
rhs_fxy = cell(nPolys_arr_fxy-1,1);

% Get RHS vector of polynomials f_{0} ... f_{d-1}
for i = 1:1:nPolys_arr_fxy-1;
    
    % Get coefficients of f(x,y)
    fxy = arr_fxy{i,1};
    
    % Get degree of f(x,y)
    m = vDeg_arr_fxy(i);
    
    % Get number of coefficients in f(x,y)
    nCoefficients_fxy = nchoosek(m+2,2);
    
    % Get the nonzero coefficients of f(x,y) from upper left triangle of
    % the coefficient matrix.
    v_fxy = GetAsVector(fxy);
    rhs_fxy{i,1} = v_fxy(1:nCoefficients_fxy);
    
end

% Build the rhs vector
rhs_fxy = cell2mat(rhs_fxy);
end

function arr_zxy = GetArray(v_zx,v_deg_arr_fxy)
% Given the vector of perturbations of f(x) given by v_zx

% Get number of polynomials in arr_fx
nPolys_fx = size(v_deg_arr_fxy,1);

% Initialise an array
arr_zxy = cell(nPolys_fx,1);

for i = 1:1:nPolys_fx
    
    % Get degree of f_{i}(x)
    m = v_deg_arr_fxy(i);
    
    % Get number of coefficients in f_{i}(x)
    nCoefficients_fxy = nchoosek(m+2,2);
    
    % Get the coefficients from the vector
    temp_vec = v_zx(1:nCoefficients_fxy);
    
    % Remove the m+1 coefficients
    v_zx(1:v_deg_arr_fxy(i)+1) = [];
    
    % Get number of zeros in f_{i}(x,y) to be added to form a matrix
    try
        nZeros = nchoosek(m+1,2);
    catch
        nZeros = 0;
    end
    
    % Add zeros
    temp_vec = ...
        [
        temp_vec;
        zeros(nZeros,1)
        ];
    
    % Get as matrix
    mat = GetAsMatrix(temp_vec,m,m);
    
    arr_zxy{i} = mat;
    
    
end


end

function arr_hxy = Get_hxy(arr_pxy,unique_vMult)

nPolys_arr_pxy = size(arr_pxy,1);

% Initialise a count
count = 1;
for i = 1:1:nPolys_arr_pxy
    
    if i == 1
        nReps = unique_vMult(i);
    else
        nReps = (unique_vMult(i) - unique_vMult(i-1));
    end
    
    for j = 1:1:nReps
        arr_hxy{count,1} = arr_pxy{i};
        count = count + 1;
    end
    
end
end