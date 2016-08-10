function [arr_hxy] = Deconvolve_Bivariate_Batch_Constrained(arr_fxy,vDegt_fxy)
%
%
%
% % Inputs.
% 
% arr_fxy : Array of polynomials f(x,y) in Bernstein form.
%
% vDegt_fxy : Vector of degrees of the polynomials f_{i}(x,y).


vDeg_f = vDegt_fxy;
vDeg_h = diff(vDeg_f);
vDeg_wxy = diff([vDeg_h; 0]);
vMult = find(vDeg_wxy~=0);


% Get number of polynomials in arr_fxy
nPolys_fxy = size(arr_fxy,1);

% %
% %
% Build the LHS Matrix
LHS_Matrix = BuildDTQ(arr_fxy,vDegt_fxy);

% %
% %
% Build the RHS Vector

rhs_array = cell(nPolys_fxy-1,1);

% Get RHS vector of polynomials f_{0} ... f_{d-1}
for i = 1:1:nPolys_fxy-1;
    
    % Get coefficients of f(x,y)
    fxy = arr_fxy{i,1};
    
    % Get degree of f(x,y)
    m = vDegt_fxy(i);
    
    % Get number of coefficients in f(x,y)
    nCoefficients_fxy = nchoosek(m+2,2);
    
    % Get the nonzero coefficients of f(x,y) from upper left triangle of 
    % the coefficient matrix.
    v_fxy = GetAsVector(fxy);
    rhs_array{i,1} = v_fxy(1:nCoefficients_fxy);
            
end

% Build the rhs vector
rhs_vec = cell2mat(rhs_array);


x_ls = SolveAx_b(LHS_Matrix,rhs_vec);

unique_vMult = unique(vMult);

arr_px = cell(length(unique_vMult),1);

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
    
    arr_px{i,1} = GetAsMatrix(vec_px,deg_px,deg_px);
    
    
    
end

%
%
%
%
nEntries_px = length(arr_px);

count = 1;
for i = 1:1:nEntries_px
        
    if i == 1
        nReps = unique_vMult(i);
    else
        nReps = (unique_vMult(i) - unique_vMult(i-1));
    end
        
    for j = 1:1:nReps
        arr_hxy{count,1} = arr_px{i};
        count = count + 1; 
    end
    
end

end


function LHS_Matrix = BuildDTQ(arr_fxy,vDegt_fxy)

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
    
    Cf{i} = [];
    
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
        Tf{j} = BuildT1(fxy,deg_fx,deg_hx);
        
        % Build the diagonal matrix D
        D{j} = BuildD(deg_fx,deg_hx);
        
        % Stack beneath all other T_f
        Cf{i} = [Cf{i} ; D{j}*Tf{j}];
        
    end
    
    Q{i} = BuildQ1(deg_hx);
    
    DTQ{i} = Cf{i} * Q{i};
end

LHS_Matrix = blkdiag(DTQ{:});

end