function [] = plotMaxMinRowNorms(vMaxRowNormR1, vMinRowNormR1, limits_k, limits_t)
%
% % Inputs
%
% vMaxRowNormR1 : (Vector)
% 
% vMinRowNormR1 : (Vector)
%
% limits_k : (Int Int)
% 
% limits_t : (Int Int)

% Get upper and lower limit of k
lowerLimit_k = limits_k(1);
upperLimit_k = limits_k(2);

% Get upper and lower limit of t
lowerLimit_t = limits_t(1);
upperLimit_t = limits_t(2);


vec_x = lowerLimit_k:1:upperLimit_k;

vMetric = vMinRowNormR1 ./ vMaxRowNormR1;


figure_name = '';
figure('name',figure_name)
hold on
plot(vec_x, log10(vMetric));

vline(lowerLimit_t);
vline(upperLimit_t);

hold off



end
