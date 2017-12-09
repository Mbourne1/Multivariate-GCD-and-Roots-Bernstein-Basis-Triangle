function [] = Experiment_DeconvolutionMethod(ex_num)

close all; clc; 
global SETTINGS
SETTINGS.PLOT_GRAPHS = false;


% % Good Examples
%
% 1 : Small example, Great results for batch methods (1e-8) 
%
% 2 : Good example, shows good effect of preproc (1e-6)
%
% 3 : 
%
% 4 : 
%
% % Bad Examples
% 
% 5
%
%
%

emin = 1e-6;

bool_preproc = false;            
o_deconvolution(ex_num, emin, bool_preproc)

bool_preproc = true;            
o_deconvolution(ex_num, emin, bool_preproc)



end