function [arr_hxy] = Deconvolve_Bivariate_Batch_Constrained_With_STLN(arr_fxy,vDegt_fxy)
%
%
%
% % Inputs.
%
% arr_fxy : Array of polynomials f(x,y) in Bernstein form.
%
% vDegt_fxy : Vector of degrees of the polynomials f_{i}(x,y).


vDeg_arr_fxy = vDegt_fxy;
vDeg_h = diff(vDeg_arr_fxy);
vDeg_wxy = diff([vDeg_h; 0]);
vMult = find(vDeg_wxy~=0);


% Get the number of polynomials in the array of polynomials f_{i}(x)
nPolys_arr_fxy = size(arr_fxy,1);

% Get the number of polynomials in the array of polynomials h_{i}(x)
nPolys_arr_hxy = nPolys_arr_fxy - 1;

% % Get the degree of each of the polynomials f_{i}(x,y)

% Initialise vector
vDeg_arr_fxy = zeros(nPolys_arr_fxy,1);

% For each polynomial f_{i} get the degree
for i = 1:1:nPolys_arr_fxy
    vDeg_arr_fxy(i) = GetDegree(arr_fxy{i});
end

% For each polynomial h_{i}
vDeg_arr_hxy = vDeg_arr_fxy(1:end-1) - vDeg_arr_fxy(2:end);

%
% y - Preprocess
% n - Dont preprocess
%

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
    arr_fww{i} = GetWithThetas(arr_fxy{i},vDeg_arr_fxy(i),th1,th2);
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
v_pww = x_ls;
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
nPolys_arr_pxy = size(arr_pww,1);
vDeg_arr_pxy = zeros(nPolys_arr_pxy,1);
for i = 1:1: nPolys_arr_pxy
    vDeg_arr_pxy(i,1) = GetDegree(arr_pww{i});
end

arr_hww = Get_hxy(arr_pww,unique_vMult);



% %
% %
% Build the array of polynomails z(\omega) which are the structured
% perturbations of the array of polynomials f(x).
arr_zww = cell(nPolys_arr_fxy,1);
for i = 1:1:nPolys_arr_fxy
    
    % Get degree of polynomial f_{i}(x,y)
    m = vDeg_arr_fxy(i);
    
    arr_zww{i} = zeros(m+1,m+1);
    
end

% Build vector z(\omega) consisting of all vectors in z_{i}(x)
v_zww = Get_v_zxy(arr_zww);

% %
% %
% Build the matrix P

% Get number of coefficients across all polynomials f_{i}(x,y)
v_nCoefficients_arr_fxy = zeros(nPolys_arr_fxy,1);

for i = 1:1:nPolys_arr_fxy
    m = vDeg_arr_fxy(i);
    v_nCoefficients_arr_fxy(i,1) = nchoosek(m+2,2);
end

% Get number of coefficients across all polynomials h_{i}(x,y)
v_nCoefficients_arr_hxy = zeros(nPolys_arr_hxy,1);

for i = 1:1:nPolys_arr_hxy
    n = vDeg_arr_hxy(i);
    v_nCoefficients_arr_hxy(i,1) = nchoosek(n+2,2);
end

% Get number of coefficients across all polynomials p_{i}(x,y)
v_nCoefficients_arr_pxy = zeros(nPolys_arr_pxy,1);
for i = 1:1:nPolys_arr_pxy
    n = vDeg_arr_pxy(i);
    v_nCoefficients_arr_pxy(i,1) = nchoosek(n+2,2);
end

% Get number of coefficients in all polynomials f_{i}(x,y)
nCoefficients_fxy = sum(v_nCoefficients_arr_fxy);

nCoefficients_hxy = sum(v_nCoefficients_arr_hxy);

nCoefficients_pxy = sum(v_nCoefficients_arr_pxy);

% Get sum of all coefficients except the final f_{i}(x,y)
nCoefficients_rhs = sum(v_nCoefficients_arr_fxy(1:end-1));

% Finally build the matrix P
P = [...
    eye(nCoefficients_rhs) ...
    zeros(nCoefficients_rhs,v_nCoefficients_arr_fxy(end))
    ];


% DY_hQ
DY_hQ = BuildDYQ(arr_hww,vDeg_arr_fxy);

% Set iteration number
ite = 1;

% Build the matrix F
F = eye(nCoefficients_pxy + nCoefficients_fxy);

%
%
%
% Build the matrix G

% Build component H_h of G
H_h = DT_fwwQ;

% Build component H_z of G
H_z = DY_hQ - P;

G = [H_h H_z];

% Compute the first residual
res_vec = rhs_fww + (P*v_zww) - (DT_fwwQ * v_pww);

% Update Matrix P*z
Pz = P*v_zww;

% Perform test
% for i = 1 : 1 : nPolys_arr_fx
%     vec_fw = [RHS_vec_fww; arr_fx{i}];
% end
%test1 = DY_hQ * vec_fw;
%test2 = DT_fwQ * v_pw;
%test1./test2

condition(ite) = norm(res_vec)./norm(rhs_fww + Pz);

start_point = ...
    [
    v_pww;
    v_zww;
    ];

%Get the iterated value
yy = start_point;

% Get
s = -(yy - start_point);


while (condition(ite) > SETTINGS.MAX_ERROR_DECONVOLUTIONS)  && ...
        (ite < SETTINGS.MAX_ITERATIONS_DECONVOLUTIONS)
    
    % Use the QR decomposition to solve the LSE problem and then
    % update the solution.
    % min |Fy-s| subject to Gy=t
    y = LSE(F,s,G,res_vec);
    
    yy = yy + y;
    
    % output y gives delta p and delta z
    delta_pww = y(1:nCoefficients_pxy);
    delta_zww = y(nCoefficients_pxy + 1:end);
    
    % Add structured perturbations to vector p(\omega)
    v_pww = v_pww + delta_pww;
    
    % Add structured perturbations to vector z(\omega)
    v_zww = v_zww + delta_zww;
    
    % Get the updated array of polynomials p_{i}(\omega)
    arr_pww = GetArray(v_pww, vDeg_arr_pxy);
    
    arr_zww = GetArray(v_zww, vDeg_arr_fxy);
    
    arr_hww = Get_hxy(arr_pww,unique_vMult);
    
    s = -(yy-start_point);
    
    DY_hQ = BuildDYQ(arr_hww,vDeg_arr_fxy);
    
    % Build the matrix C(f)
    DT_fwwQ = BuildDTQ(arr_fww,vDeg_arr_fxy);
    
    % Build the matrix C(z)
    DT_zwwQ = BuildDTQ(arr_zww,vDeg_arr_fxy);
    
    % Build G
    H_z = DY_hQ - P;
    H_h = DT_fwwQ + DT_zwwQ;
    
    G = [H_h H_z];
    
    % Update the RHS vector
    rhs_fww = BuildRHS_vec(arr_fww);
    RHS_vec_Pzww = BuildRHS_vec(arr_zww);
    
    
    % Calculate residual and increment t in LSE Problem
    res_vec = ((rhs_fww + RHS_vec_Pzww) - ((DT_fwwQ + DT_zwwQ)*v_pww));
    
    % Increment iteration number
    ite = ite + 1;
    
    % Get condition number
    condition(ite) = norm(res_vec)./norm(rhs_fww + RHS_vec_Pzww);
end


% Plot condiiton

figure()
plot(log10(condition))
hold off


% Get array of polynomials h_{i}(x) from h_{i}(\omega)
arr_hxy = cell(nPolys_arr_hxy,1);
for i = 1:1:nPolys_arr_hxy
    m = vDeg_arr_hxy(i);
    arr_hxy{i} = GetWithoutThetas(arr_hww{i},m,th1,th2);
end



end


function LHS_Matrix = BuildDTQ(arr_fxy,vDegt_fxy)

vDeg_f = vDegt_fxy;
vDeg_h = diff(vDeg_f);
vDeg_wxy = diff([vDeg_h; 0]);
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


function rhs_vec = BuildRHS_vec(arr_fxy)

% Get number of polynomials in array of f_{i}(x,y)
nPolys_arr_fxy = size(arr_fxy,1);

% Get the degree of the polynomials
vDeg_arr_fxy = zeros(nPolys_arr_fxy,1);
for i = 1:1:nPolys_arr_fxy
    vDeg_arr_fxy(i) = GetDegree(arr_fxy{i});
end


% Get RHS array
rhs_array = cell(nPolys_arr_fxy-1,1);

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
    rhs_array{i,1} = v_fxy(1:nCoefficients_fxy);
    
end

% Build the rhs vector
rhs_vec = cell2mat(rhs_array);
end


function v_zxy =  Get_v_zxy(arr_zxy)

% Get the number of polynomials in array
nPolys_arr_zxy = size(arr_zxy,1);

% Get the degree of each of the polynomials z(x,y)
vDeg_arr_zxy = zeros(nPolys_arr_zxy);

for i = 1:1:nPolys_arr_zxy
    vDeg_arr_zxy(i) = GetDegree(arr_zxy{i});
end

% For each polynomial in the array, get as a vector
for i = 1:1:nPolys_arr_zxy
    
    % Get as vector
    v_zxy = GetAsVector(arr_zxy{i});
    
    % Get number of non-zeros
    nZeros_zxy = nchoosek(vDeg_arr_zxy(i)+2,2);
    
    % Remove zeros from vector
    v_zxy = v_zxy(1:nZeros_zxy);
    
    arr_zxy{i} = v_zxy;
end

v_zxy = cell2mat(arr_zxy);

end

function arr_hxy = Get_hxy(arr_pxy,vUniqueMult)

% Get number of entries in the array of polynomials p_{i}(x)
nPolys_arr_px = size(arr_pxy,1);

% initialise count
count = 1;


for i = 1:1:nPolys_arr_px
    
    if i == 1
        nReps = vUniqueMult(i);
    else
        nReps = (vUniqueMult(i) - vUniqueMult(i-1));
    end
    
    for j = 1:1:nReps
        arr_hxy{count,1} = arr_pxy{i};
        count = count + 1;
    end
    
end
end


function C_fxy = BuildDYQ(arr_hxy,vDeg_arr_fxy)
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
nPolys_arr_hxy = size(arr_hxy,1);

% Get degree of each polynomial f_{i}(x,y)
vDeg_arr_hxy = zeros(nPolys_arr_hxy,1);
for i = 1 : 1 : nPolys_arr_hxy
    vDeg_arr_hxy(i) = GetDegree(arr_hxy{i});
end


% Initialise a cell array
arr_DT1Q1 = cell(nPolys_arr_hxy,1);

% For each of the polynomials excluding the first f_{1},...,f_{d}
for i = 1:1:nPolys_arr_hxy
    
    
    % Get the degree of h{i}(x,y)
    m = vDeg_arr_hxy(i);
    
    % Get the degree of f{i}(x,y)
    n = vDeg_arr_fxy(i+1);
    
    % Temporarily call the ith entry f(x,y)
    hxy = arr_hxy{i};
    
    % Build the matrix T_{n-m}(f(x,y))
    T1 = BuildT1(hxy,m,n);
    
    D = BuildD(m,n);
    Q1 = BuildQ1(n);
    
    arr_DT1Q1{i} = D*T1*Q1;
end

% Include a section of zeros

C_fxy = blkdiag(arr_DT1Q1{:});
nRows = size(C_fxy,1);

nCols = nchoosek(vDeg_arr_fxy(1)+2,2);
zero_section = zeros(nRows,nCols);

C_fxy = [zero_section C_fxy];
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
    v_zx(1:nCoefficients_fxy) = [];
    
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