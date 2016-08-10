function [] = SetGlobalVariables(mean_method,bool_alpha_theta)

global SETTINGS

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

SETTINGS.THRESHOLD = 2;
SETTINGS.THRESHOLD_RANK = 1e-10;

%
% 'y'
% 'n'
%
SETTINGS.PLOT_GRAPHS = 'y';

SETTINGS.SEED = 1024;

%
% 'Separate'
% 'Batch'
% 'Batch Constrained'
%
%
SETTINGS.DECONVOLUTION_METHOD = 'Batch Constrained';
end
