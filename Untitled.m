

% Example corresponds to example in report on gcd computation
ex_num = '12';
el = 0;
eu = 0;

close all; clc; 

o_gcd_Bivariate_2Polys(ex_num, el, eu, 'None', false, 'None', 'None', 'DTQ')
o_gcd_Bivariate_2Polys(ex_num, el, eu, 'Geometric Mean Matlab Method', true, 'None', 'None', 'DTQ')