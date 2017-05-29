function [] = SetGlobalVariables( ex_num, emin, emax,...
    mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, sylvester_build_method)
%
% % Inputs
%
% ex_num : (String) Example Number
%
% emin : (Float)
%
% emax : (Float)
%
% mean_method : (String)
%
% bool_alpha_theta : (Boolean)
%
% low_rank_approx_method : (String)
%
% apf_method : (String)
%
% sylvester_build_method : (String)

global SETTINGS

%----------------------

% Plot Graphs
%   * true
%   * false
%
SETTINGS.PLOT_GRAPHS = true;



SETTINGS.SEED = 1024;

% SYLVESTER BUILD METHOD
%   T
%   DTQ
%   DT
%   TQ
SETTINGS.SYLVESTER_BUILD_METHOD = sylvester_build_method;

% Set Noise
SETTINGS.EMIN = emin;
SETTINGS.EMAX = emax;

% Set example number
SETTINGS.EX_NUM = ex_num;

%--------------------------------------------------------------------------
%
%       SETTINGS - PREPROCESSING
%
%


%
% 'Geometric Mean Matlab Method'
% 'Geometric Mean My Method'
% 'None'
%
SETTINGS.MEAN_METHOD = mean_method;

%
%   * true
%   * false
%
SETTINGS.BOOL_ALPHA_THETA = bool_alpha_theta;



% -------------------------------------------------------------------------
%
%   SETTINGS - DEGREE COMPUTATION
%
%
SETTINGS.THRESHOLD = 2;
SETTINGS.THRESHOLD_RANK = 1e-3;

% Metric used to compute the degree of the GCD
% R1 Row Norms
% R1 Row Diagonals
% Singular Values
% Residuals
SETTINGS.RANK_REVEALING_METRIC = 'Singular Values';

% ------------------------------------------------------------------------
%
%   SETTINGS - LOW RANK APPROXIMATION
%
%


% 'Standard STLN'
% 'None'
%
SETTINGS.LOW_RANK_APPROX_METHOD = low_rank_approx_method;

SETTINGS.STLN_MAX_ERROR = 1e-10;
SETTINGS.STLN_MAX_ITERATIONS = 50;

%--------------------------------------------------------------------------
%
%   SETTINGS - APPROXIMATE FACTORISATION
%
%
SETTINGS.APF_METHOD = apf_method;


%-------------------------------------------------------------------------
%
%   SETTINGS - DECONVOLUTION
%
%
%

% Deconvolution Settings for root finding method

% Method used to deconvolve polynomials f_{i}(x,y) to compute h(x,y)
%
% 'Separate'
% 'Batch Without STLN'
% 'Batch With STLN'
% 'Batch Constrained Without STLN'
% 'Batch Constrained With STLN'
%
SETTINGS.HXY_DECONVOLUTION_METHOD = 'Batch Constrained Without STLN';


% Method used to deconvolve polynomials h_{i}(x,y) to compute w_{i}(x,y)
%
% 'Separate'
% 'Batch Without STLN'
% 'Batch With STLN'
%
SETTINGS.WXY_DECONVOLUTION_METHOD = 'Batch Without STLN';


SETTINGS.PREPROC_DECONVOLUTIONS = true;
end
