function [] = Experiment2Preprocessing_2Polys(ex_num)
%
%
% >> Experimtne2Preprocessing_2Polys('19')


addpath(genpath(pwd));

close all; 
clc;

% Good Examples



% Constants

el = 1e-14;
eu = 1e-14;

low_rank_approx_method = 'None';
apf_method = 'None';
rank_revealing_metric = 'Minimum Singular Values';

sylvester_format = 'DTQ';



mean_method = 'None';
bool_preproc = false;

o_gcd_Bivariate_2Polys(ex_num, el, eu, mean_method, bool_preproc, low_rank_approx_method, apf_method, sylvester_format, rank_revealing_metric)

mean_method = 'Geometric Mean Matlab Method';
bool_preproc = true;

o_gcd_Bivariate_2Polys(ex_num, el, eu, mean_method, bool_preproc, low_rank_approx_method, apf_method, sylvester_format, rank_revealing_metric)


end
