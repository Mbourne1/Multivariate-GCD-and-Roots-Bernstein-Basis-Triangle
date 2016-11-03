function [] = SetGlobalVariables(mean_method,bool_alpha_theta,low_rank_approx_method,sylvester_type)

global SETTINGS

%----------------------

% Plot Graphs
% 'y'
% 'n'
%
SETTINGS.PLOT_GRAPHS = 'y';

SETTINGS.SEED = 1024;

% T, DTQ, DT, TQ
SETTINGS.SYLVESTER_MATRIX_TYPE = sylvester_type;

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
% 'y'
% 'n'
%
SETTINGS.BOOL_ALPHA_THETA = bool_alpha_theta;



% -------------------------------------------------------------------------
%
%   SETTINGS - DEGREE COMPUTATION
%
%
SETTINGS.THRESHOLD = 2;
SETTINGS.THRESHOLD_RANK = 1e-10;


% ------------------------------------------------------------------------
%
%   SETTINGS - LOW RANK APPROXIMATION
%
%


% 'Standard STLN'
% 'None'
%
SETTINGS.LOW_RANK_APPROX_METHOD = low_rank_approx_method;

SETTINGS.STLN_MAX_ERROR = 1e-13;
SETTINGS.STLN_MAX_ITERATIONS = 100;

%--------------------------------------------------------------------------


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


SETTINGS.PREPROC_DECONVOLUTIONS = 'y';
end
