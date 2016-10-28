function [] = o_roots(ex_num,emin,mean_method,bool_preproc,low_rank_approx_method)
% O_ROOTS
% O_ROOTS(ex_num,emin,emax,mean_method,bool_alpha_theta,low_rank_approx_method)
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
% >> o_roots('1', 1e-10, 'Geometric Mean Matlab Method', 'y','Standard STLN')

% Set variables
global SETTINGS
SETTINGS.PLOT_GRAPHS = 'n';

% % 
% Add subfolders
addpath(...
    'Deconvolution',...
    'Examples',...
    'Examples/Examples_Roots',...
    'Formatting',...
    'Get Cofactors',...
    'Get GCD Coefficients',...
    'Get GCD Degree',...    
    'Preprocessing',...
    'Root Finding Methods',...
    'Root Finding Methods/My Method'...
);

SetGlobalVariables(mean_method,bool_preproc,low_rank_approx_method);

% Set the problem type, used in logging to the outputs file.
problem_type = 'Roots Bivariate Bernstein Triangle';

% %
% % 
% Get Example

% Given the example number get the polynomial coefficients of f(x,y) and
% the total degree 'm' of f(x,y).
[fxy_exact,m] = Examples_Roots(ex_num);

% Add noise to the coefficients of polynomial f(x,y)
[fxy,~] = Noise2(fxy_exact,emin);

o_roots_mymethod(fxy,m);


end