function [] = o_roots(ex_num, emin, emax, mean_method, bool_alpha_theta,...
    low_rank_approx_method, apf_method, sylvester_matrix_type)
% O_ROOTS(ex_num,emin,emax,mean_method,bool_alpha_theta,low_rank_approx_method, apf_method, sylvester_matrix_type)
%
% Given an example number and set of parameters, obtain the roots of the
% example polynomial, where the polynomial is in the power basis
%
% %                             Inputs
%
% ex_num - Example Number
%
% emin : Lower noise level
%
% emax : Upper noise level
%
% mean_method :
%       'Geometric Mean Matlab Method'
%       'None'
%
% bool_alpha_theta ('y'/'n')
%       'y' - Include Preprocessing
%       'n' - Exclude Preprocessing
%
% low_rank_approx_method ('y'/'n')
%       'Standard SNTLN' : Include SNTLN
%       'None'           : Exclude SNTLN
%
% % Examples
% >> o_roots('1', 1e-12, 1e-10, 'Geometric Mean Matlab Method', 'y','Standard STLN','Standard APF','DTQ')

% Set variables
global SETTINGS
SETTINGS.PLOT_GRAPHS = 'n';

% %
% Add subfolders
addpath(...
    'APF',...
    'Bernstein Functions',...
    'Build Matrices',...
    'Build Sylvester Matrix',...
    'Deconvolution',...
    'Examples',...
    'Examples/Examples_Roots',...
    'Formatting',...
    'GCD Methods',...
    'Get Cofactors',...
    'Get GCD Coefficients',...
    'Get GCD Degree',...
    'Low Rank Approximation',...
    'Preprocessing',...
    'Root Finding Methods',...
    'Root Finding Methods/My Method'...
    );

% Set the problem type, used in logging to the outputs file.
problem_type = 'Roots Bivariate Bernstein Triangle';

% Set the deconvolution method for the batch of deconvolutions in the
% factorisation algorithm
SETTINGS.DECONVOLUTION_METHOD = 'Batch';

SetGlobalVariables(problem_type, ex_num, emin, emax, mean_method, ...
    bool_alpha_theta, low_rank_approx_method, apf_method, sylvester_matrix_type);


% %
% %
% Get Example

% Given the example number get the polynomial coefficients of f(x,y) and
% the total degree 'm' of f(x,y).
[fxy_exact,m] = Examples_Roots(ex_num);

% Add noise to the coefficients of polynomial f(x,y)
[fxy,~] = AddNoiseToPoly(fxy_exact,emin);


% %

root_finding_method = '3 Poly GCD';

switch root_finding_method
    case '2 Poly GCD'
        
        o_roots_mymethod(fxy,m);
        
    case '3 Poly GCD'
        
        o_roots_mymethod_newmethod(fxy,m);
        
    otherwise
        error([mfilename ': Error']);
end
end