function [] = Test_Deconvolution(ex_num,emin)
% Perform tests on deconvolution methods
%
%
% % Inputs
%
% ex_num : Example number (String)
%
% emin : Noise level
%
%
% % Outputs.
%
%  
%
% Example
%
% >> Test_Deconvolution('1',1e-10);



syms x y;

global SETTINGS
SETTINGS.SEED = 1024;


switch ex_num
    case '1' 
        
        % (x + y + 0.5)^7 * (x+y+1.2)^12
        
        factor(1,1) = (x + y + 0.5);
        factor(2,1) = (x + y + 1.2);
        
        % Set the multiplicity of each factor
        vMult = [7; 12];
        

    case '2'
        factor(1,1) = (x + y + 0.5);
        factor(2,1) = (x + y + 1.2);
        
        % Set the multiplicity of each factor
        vMult = [10; 15];
        
end


% Get the highest power of any factor
highest_pwr = max(vMult);


% %
% %
% Generate the polynomials f_{0},...,f_{m}(x)

arr_sym_fxy = cell(highest_pwr+1,1);
vDegt_arr_fxy = zeros(highest_pwr+1,1);
for i = 0:1:highest_pwr
    
    % Get the multiplicities of the roots of f_{i+1}
    mults = ((vMult - i) + abs(vMult-i)) ./ 2;
    
    % Get the symbolic polynomial f_{i+1}
    arr_sym_fxy{i+1} = prod(factor.^(mults));
    
    % Get the degree of polynomial f_{i+1}(x)
    vDegt_arr_fxy(i+1) = double(feval(symengine, 'degree', (arr_sym_fxy{i+1})));
end

% Get the degree structure of the polynomials h_{i} where h_{i} =
% f_{i-1}/f_{i}.
vDegt_arr_hxy = -1.*diff(vDegt_arr_fxy);

% Get the degree structure of the polynomials w_{i}
vDegt_arr_wxy = diff([vDegt_arr_hxy; 0]);

% Get the multiplicities of the roots
vMultiplicities = find(vDegt_arr_wxy~=0);

% Get the sequence of polynomials h_{i}(x) in symbolic form
arr_sym_h = cell(length(arr_sym_fxy)-1,1);
for i = 1:1:length(arr_sym_fxy)-1
    arr_sym_h{i} = arr_sym_fxy{i} / arr_sym_fxy{i+1};
end


% %
% %
% Get coefficients vectors of f_{i}(x) and h_{i}(x)
nPolys_arr_fxy = length(arr_sym_fxy);
nPolys_arr_hxy = length(arr_sym_fxy) - 1;

arr_fxy_pwr = cell(nPolys_arr_fxy,1);
arr_fxy_brn = cell(nPolys_arr_fxy,1);

arr_hxy_pwr = cell(nPolys_arr_hxy,1);
arr_hxy_brn = cell(nPolys_arr_hxy,1);

for i = 1:1:nPolys_arr_fxy
    
        
        arr_fxy_pwr{i,1} = double(rot90(coeffs(arr_sym_fxy{i},[x,y],'All'),2));
        
        arr_fxy_brn{i,1} = PowerToBernstein(arr_fxy_pwr{i,1},vDegt_arr_fxy(i));
        
    if i <= nPolys_arr_hxy    
        arr_hxy_pwr{i,1} = double(rot90(coeffs(arr_sym_h{i},[x,y],'All'),2));
        
        arr_hxy_brn{i,1} = PowerToBernstein(arr_hxy_pwr{i,1},vDegt_arr_hxy(i));
        
    else
        arr_fxy_pwr{i,1} = 1;
    end
    
end

% %
% %
% %
% Add noise to the coefficients of f_{i}(x)
arr_fxy_brn_noisy = cell(nPolys_arr_fxy,1);

for i = 1:1:nPolys_arr_fxy
    arr_fxy_brn_noisy{i,1} = Noise(arr_fxy_brn{i},emin);
end

% -------------------------------------------------------------------------
% %
% %
% %
% Testing deconvolution batch method with constraints

LineBreakLarge()
fprintf('Deconvolve with the constraint that all h_{i} for i between m_{j} and m_{j+1} are equal \n')
fprintf('Staggered Staircase Method \n\n')

arr_hx_test_1_brn = Deconvolve_Bivariate_Batch_Constrained(arr_fxy_brn_noisy,vDegt_arr_fxy);

for i = 1:1:length(arr_hx_test_1_brn)
    %display(arr_hx_test_1{i} );
end

% Compare each computed h{i} with actual h_{i}
for i = 1:1:length(arr_hx_test_1_brn)
    
    exact = arr_hxy_brn{i}./arr_hxy_brn{i}(1,1);
    comp = arr_hx_test_1_brn{i}./arr_hx_test_1_brn{i}(1,1);
    
    err_measure = norm(exact - comp) ./ norm(exact);
    fprintf([mfilename ': ' sprintf('%i : Error : %2.4e \n', i,err_measure)]);
end

% %
% %
% %
% Testing deconvolution batch method with constraints and low rank approx

LineBreakLarge()
fprintf('Deconvolve with the constraint that all h_{i} for i between m_{j} and m_{j+1} are equal \n')
fprintf('With STLN \n')
fprintf('Staggered Staircase Method \n\n')

arr_hxy_brn_test_3 = Deconvolve_Bivariate_Batch(arr_fxy_brn_noisy,vDegt_arr_fxy);

for i = 1:1:length(arr_hxy_brn_test_3)
    %display(arr_hx_test_1{i} );
end

% Compare each computed h{i} with actual h_{i}
for i = 1:1:length(arr_hxy_brn_test_3)
    
    exact = arr_hxy_brn{i}./arr_hxy_brn{i}(1,1);
    comp = arr_hxy_brn_test_3{i}./ arr_hxy_brn_test_3{i}(1,1);
    
    err_measure = norm(exact - comp) ./ norm(exact);
    fprintf([mfilename ': ' sprintf('%i : Error : %2.4e \n', i,err_measure)]);
end




%--------------------------------------------------------------------------
% %
% %
% %
% Testing deconvolution
LineBreakLarge();
fprintf('Deconvolve - Each deconvolution is independent \n');
fprintf('Deconvolution Separate \n');


arr_hxy_brn_test_4 = Deconvolve_Bivariate_Separate(arr_fxy_brn_noisy);

for i = 1:1:length(arr_hxy_brn_test_4)
    %display(arr_hx_test_3{i});
end
% Compare each computed h{i} with actual h_{i}
for i = 1:1:length(arr_hxy_brn_test_4)
    
    exact = arr_hxy_brn{i}./arr_hxy_brn{i}(1,1);
    comp = arr_hxy_brn_test_4{i}./arr_hxy_brn_test_4{i}(1,1);
    
    err_measure = norm(exact - comp) ./ norm(exact);
    fprintf([mfilename ': ' sprintf('%i : Error : %2.4e \n', i,err_measure)]);
end


end

