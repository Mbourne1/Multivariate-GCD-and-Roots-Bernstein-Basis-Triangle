function [] = Experiment3(ex_num, bool_preproc)
close all; clc;

% Examples
%
% 1
% 2
% 3 - bad example

close all; clc;


emin = 1e-10;
emax = 1e-8;


switch bool_preproc
    
    case true
        mean_method = 'Geometric Mean Matlab Method';
        bool_alpha_theta = true;
        
    case false
        mean_method = 'None';
        bool_alpha_theta = false;
end

low_rank_approx_method = 'None';
apf_method = 'None';
sylvester_matrix_type = 'DTQ';
rank_revealing_method = 'Minimum Singular Values';
deconvolution_method_wxy = 'Batch';

nEquations = '2';

% Experiment 1

arr_deconvolution_method_hxy = {'Separate', 'Batch',  'Batch With STLN', 'Batch Constrained', 'Batch Constrained With STLN'};
nDeconvolutionMethods = length(arr_deconvolution_method_hxy);

for i = 1 : 1 : nDeconvolutionMethods
    
    deconvolution_method_hxy = arr_deconvolution_method_hxy{i};
    
    
    o_roots_Bivariate(ex_num, emin, emax, mean_method, bool_alpha_theta, ...
        low_rank_approx_method, apf_method, sylvester_matrix_type, ...
        rank_revealing_method, deconvolution_method_hxy, deconvolution_method_wxy, nEquations)
    
end

end