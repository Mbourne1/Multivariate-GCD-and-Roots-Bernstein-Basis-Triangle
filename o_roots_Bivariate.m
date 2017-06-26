function [] = o_roots_Bivariate(ex_num, emin, emax, mean_method, ...
    bool_alpha_theta, low_rank_approx_method, apf_method, ...
    sylvester_matrix_type, rank_revealing_metric, deconvolution_method_hxy, ...
    deconvolution_method_wxy)
% O_ROOTS_BIVARIATE(ex_num,emin,emax,mean_method,bool_alpha_theta,low_rank_approx_method, apf_method, sylvester_matrix_type)
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
    deconvolution_method_hxy, deconvolution_method_wxy);


% %
% %
% Get Example

% Given the example number get the polynomial coefficients of f(x,y) and
% the total degree 'm' of f(x,y).
[fxy_exact,m] = Examples_Roots(ex_num);

% Add noise to the coefficients of polynomial f(x,y)
[fxy,~] = AddNoiseToPoly(fxy_exact, emin);


% %

root_finding_method = '3 Poly GCD';

switch root_finding_method
    case '2 Poly GCD'
        
        o_roots_mymethod(fxy, m);
        
    case '3 Poly GCD'
        
        o_roots_mymethod_newmethod(fxy, m);
        
    otherwise
        error([mfilename ': Error']);
end
end