function [arr_hxy] = Deconvolve_Bivariate_Batch_Constrained_Without_STLN(arr_fxy,vDegt_fxy)
%
%
%
% % Inputs.
%
% arr_fxy : (Array of Matrices) Array of polynomials f(x,y) in Bernstein form.
%
% vDegt_fxy : (Vector) Vector of degrees of the polynomials f_{i}(x,y).


vDeg_arr_fxy = vDegt_fxy;
vDeg_arr_hxy = diff(vDeg_arr_fxy);
vDeg_arr_wxy = diff([vDeg_arr_hxy; 0]);
vMult = find(vDeg_arr_wxy~=0);


% Get number of polynomials in arr_fxy
nPolys_arr_fxy = size(arr_fxy,1);

% Get number of polynomials in array of h_{i}(x,y)
nPolys_arr_hxy = nPolys_arr_fxy - 1;

vDeg_arr_fxy = zeros(nPolys_arr_fxy,1);
for i = 1:1:nPolys_arr_fxy
    
    vDeg_arr_fxy(i)  = GetDegree_Bivariate(arr_fxy{i});
    
end


% Preprocess
global SETTINGS
if (SETTINGS.PREPROC_DECONVOLUTIONS)
    
    [th1, th2] = GetOptimalTheta(arr_fxy);
    
else
    
    th1 = 1;
    th2 = 1;
    
end

% Get f(\omega_{1},\omega_{2}) from f(x,y)
arr_fww = cell(nPolys_arr_fxy,1);
for i = 1:1:nPolys_arr_fxy
    
    m = vDeg_arr_fxy(i);
    arr_fww{i} = GetWithThetas(arr_fxy{i}, m, th1, th2);
    
end

% %
% %
% Build the LHS Matrix
DT_fwwQ = BuildDTQ_2Polys(arr_fww, vDeg_arr_fxy);

% %
% %
% Build the RHS Vector

rhs_fww = BuildRHS_vec(arr_fww);


x_ls = SolveAx_b(DT_fwwQ, rhs_fww);

unique_vMult = unique(vMult);

arr_pww = cell(length(unique_vMult),1);

for i = 1:1:length(unique_vMult)
    
    
    mult = unique_vMult(i);
    
    % Get degree of p(x,y)
    deg_px = vDegt_fxy(mult) - vDegt_fxy(mult+1);
    
    % Get number of coefficients in p(x,y)
    nCoefficients_px = nchoosek(deg_px+2, 2);
    
    % Get coefficients of p(x,y) from x_ls
    vec_px = x_ls(1:nCoefficients_px);
    
    % Remove coefficients
    x_ls(1:nCoefficients_px) =[];
    
    
    nZeros_px = nchoosek(deg_px+1, 2);
    
    vec_px = ...
        [
        vec_px;
        zeros(nZeros_px,1)
        ];
    
    arr_pww{i,1} = GetAsMatrix(vec_px, deg_px, deg_px);
    
    
    
end

%
%
%
%
nPolys_arr_px = length(arr_pww);

count = 1;
for i = 1 : 1 : nPolys_arr_px
    
    if i == 1
        nReps = unique_vMult(i);
    else
        nReps = (unique_vMult(i) - unique_vMult(i-1));
    end
    
    for j = 1:1:nReps
        arr_hxy{count,1} = arr_pww{i};
        count = count + 1;
    end
    
end

end


function LHS_Matrix = BuildDTQ_2Polys(arr_fxy, vDegt_fxy)
%
% % Inputs
%
% arr_fxy : (Array of Matrices)
%
% vDegt_fxy : (Array of Matrices)
%
% % Outputs
%
% LHS_Matrix : (Matrix)

vDeg_f = vDegt_fxy;
vDeg_h = diff(vDeg_f);
vDeg_wxy = diff([vDeg_h; 0]);
vMult = find(vDeg_wxy~=0);
display(vMult);

% Get number of distinct polynomials h_{i}(x)
nDistinct_hx = length(vMult);

for i = 1:1:nDistinct_hx
    
    if i>1
        old_mult = vMult(i-1);
    else
        old_mult = 0;
    end
    
    new_mult = vMult(i);
    
    arr_Cf{i} = [];
    
    for j = (old_mult+1+1) : 1 : (new_mult+1)
        
        % Get coefficients of previous f(x,y)
        fxy_prev = arr_fxy{j-1};
        
        % Get the degree of previous f(x,y)
        deg_fx_prev = vDeg_f(j-1);
        
        % Get polynomial f(x,y)
        fxy = arr_fxy{j};
        
        % Get the degree of f(x,y)
        deg_fx = vDeg_f(j);
        
        % Get the degree of polynomial h_{i}
        deg_hx = deg_fx_prev - deg_fx;
        
        % Build the cauchy like matrix
        arr_Tf{j} = BuildT1(fxy,deg_fx,deg_hx);
        
        % Build the diagonal matrix D
        arr_D{j} = BuildD(deg_fx,deg_hx);
        
        % Stack beneath all other T_f
        arr_Cf{i} = [arr_Cf{i} ; arr_D{j}*arr_Tf{j}];
        
    end
    
    arr_Q{i} = BuildQ1(deg_hx);
    
    arr_DTQ{i} = arr_Cf{i} * arr_Q{i};
end

LHS_Matrix = blkdiag(arr_DTQ{:});

end

function rhs_fxy = BuildRHS_vec(arr_fxy)
%
% % Inputs
%
% arr_fxy : (Array of Matrices)
%
% % Outputs
%
% rhs_fxy : (Vector)

% Get number of polynomials
nPolys_arr_fxy = size(arr_fxy, 1);

% Get degree of the polynomials
vDeg_arr_fxy = zeros(nPolys_arr_fxy, 1);

for i = 1:1:nPolys_arr_fxy
    
    vDeg_arr_fxy(i) = GetDegree_Bivariate(arr_fxy{i});
    
end

rhs_fxy = cell(nPolys_arr_fxy-1, 1);

% Get RHS vector of polynomials f_{0} ... f_{d-1}
for i = 1:1:nPolys_arr_fxy-1
    
    % Get coefficients of f(x,y)
    fww = arr_fxy{i, 1};
    
    % Get degree of f(x,y)
    m = vDeg_arr_fxy(i);
    
    % Get number of coefficients in f(x,y)
    nCoefficients_fxy = nchoosek(m+2, 2);
    
    % Get the nonzero coefficients of f(x,y) from upper left triangle of
    % the coefficient matrix.
    v_fxy = GetAsVector(fww);
    rhs_fxy{i,1} = v_fxy(1:nCoefficients_fxy);
    
end

% Build the rhs vector
rhs_fxy = cell2mat(rhs_fxy);

end