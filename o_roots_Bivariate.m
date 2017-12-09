function [] = o_roots_Bivariate(ex_num, emin, emax, mean_method, ...
    bool_alpha_theta, low_rank_approx_method, apf_method, ...
    sylvester_matrix_type, rank_revealing_metric, deconvolution_method_hxy, ...
    deconvolution_method_wxy, nEquations)
% O_ROOTS_BIVARIATE(ex_num,emin,emax,mean_method,bool_alpha_theta,low_rank_approx_method, apf_method, sylvester_matrix_type, rank_revealing_method, deconvolution_method_hxy, deconvolution_method_wxy)
%
% Given an example number and set of parameters, obtain the roots of the
% example polynomial, where the polynomial is in the power basis
%
% % Inputs
%
% ex_num : (String) Example Number
%
% emin : (Float) Lower noise level
%
% emax : (Float) Upper noise level
%
% mean_method : (String)
%   * 'Geometric Mean Matlab Method'
%   * 'None'
%
% bool_alpha_theta : (Boolean)
%   * true :  Include Preprocessing
%   * false :  Exclude Preprocessing
%
% low_rank_approx_method : (String)
%   * 'Standard SNTLN' : Include SNTLN
%   * 'None'           : Exclude SNTLN
%
% sylvester_build_method : (String)
%   * 'T'
%   * 'DT'
%   * 'TQ'
%   * 'DTQ'
%   * 'DTQ Denominator Removed'
%    
%
%
%
%   
%
% % Examples
% >> o_roots_Bivariate('1', 1e-12, 1e-10, 'None', false, 'None', 'None', 'DTQ', 'Minimum Singular Values', 'Batch', 'Batch')
% >> o_roots_Bivariate('1', 1e-12, 1e-10, 'Geometric Mean Matlab Method', true, 'Standard STLN', 'Standard APF', 'DTQ', 'Minimum Singular Values', 'Batch', 'Batch')



% Restore defaults and add subfolders
restoredefaultpath
addpath(genpath(pwd));


SetGlobalVariables_Roots(ex_num, emin, emax, mean_method, ...
    bool_alpha_theta, low_rank_approx_method, apf_method, ...
    sylvester_matrix_type, rank_revealing_metric, ...
    deconvolution_method_hxy, deconvolution_method_wxy, nEquations);


% %
% %
% Get Example

% Given the example number get the polynomial coefficients of f(x,y) and
% the total degree 'm' of f(x,y).
[fxy_exact, arr_fxy_exact, arr_hxy_exact, arr_wxy_exact, m] = Examples_Roots(ex_num);

% Add noise to the coefficients of polynomial f(x,y)
[fxy,~] = AddNoiseToPoly(fxy_exact, emin);


% Get the array of polynomials f_{i}(x,y), h_{i}(x,y), w_{i}(x,y)        
[arr_fxy, arr_hxy, arr_wxy, ~] = o_roots_mymethod_newmethod(fxy, m);
        

 

% Compute the distance between the array or exact polynomials f_{i}(x,y)
% and the computed polynomials f_{i}(x,y)
vErrors_arr_fxy = CompareArrays(arr_fxy_exact, arr_fxy);

% Compute the distance between the array of exact polynomials h_{i}(x,y)
% and the computed polynomials h_{i}(x,y)
vErrors_arr_hxy = CompareArrays(arr_hxy_exact, arr_hxy);

% Compute the distance between the array of exact polynomials w_{i}(x,y)
% and the computed polynomials w_{i}(x,y)
vErrors_arr_wxy = CompareArrays(arr_wxy_exact, arr_wxy);


% Plot Errors
PlotErrors({vErrors_arr_fxy, vErrors_arr_hxy, vErrors_arr_wxy},...
    {'$\epsilon f_{i}(x,y)$', ...
    '$\epsilon h_{i}(x,y)$', ...
    '$\epsilon w_{i}(x,y)$'});

% Print Errors
LineBreakLarge()
fprintf('Average Error f_{i}(x,y) : %e \n', mean(vErrors_arr_fxy))
fprintf('Average Error h_{i}(x,y) : %e \n', mean(vErrors_arr_hxy))
fprintf('Average Error w_{i}(x,y) : %e \n', mean(vErrors_arr_wxy))
LineBreakLarge()


end

function [vDistance] = CompareArrays(arr_fxy_exact, arr_fxy)


% Get length of array
nPolysArray = length(arr_fxy_exact);

% Store distance
vDistance = zeros(nPolysArray,1);

for i = 1 : 1 : nPolysArray
    
    fxy_exact = arr_fxy_exact{i};
    fxy_comp = arr_fxy{i};

    vDistance(i) = GetMatrixDistance(fxy_exact, fxy_comp);
end



end


function PlotErrors(arrErrors, arrStrings)
%
% % Input
%
% arrErrors : (Array of Vectors) Each cell of the array contains a vector
% of errors
%
%
% arrStrings : (Array of Strings)
figure_name = 'Errors';
figure('Name',figure_name)

hold on
for i = 1 : 1 : 3
   
    vErrors = arrErrors{i};
    nErrors = length(vErrors);
    x_vec = 1 : 1 : nErrors;
    strName = arrStrings{i};
    plot(x_vec, log10(vErrors), '-o', 'DisplayName', strName, 'LineWidth', 2);
    
end

l = legend(gca,'show');
set(l,'Interpreter', 'latex');

grid on
box on

hold off

end
