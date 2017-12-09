function [] = Experiment1SylvesterVariants_3Polys(ex_num)

addpath(genpath(pwd));

%close all;
%clc;

% Good Examples
% ex_num = '14'; el = 1e-7; eu = 1e-7;
% ex_num = '14'; el = 1e-6; eu = 1e-6;

% Constants
el = 1e-6;
eu = 1e-5;

mean_method = 'None';
bool_preproc = false;

%mean_method = 'Geometric Mean Matlab Method';
%bool_preproc = true;

low_rank_approx_method = 'None';
apf_method = 'None';
rank_revealing_metric = 'Minimum Singular Values';


% %

arrSylvesterFormat = {'T', 'DT', 'TQ', 'DTQ'};
nEquations = '3';

for i = 1 : 1 : length(arrSylvesterFormat)
    
    
    sylvester_format = arrSylvesterFormat{i};
    o_gcd_Bivariate_3Polys(ex_num, el, eu, mean_method, bool_preproc, ...
        low_rank_approx_method, apf_method, sylvester_format, ...
        rank_revealing_metric, nEquations)
    
end

end