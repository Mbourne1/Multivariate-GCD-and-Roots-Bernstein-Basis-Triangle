function [] = plotMaxMinRowNorms(vMaxRowNormR1, vMinRowNormR1, limits_t)
%
% % Inputs
%
% vMaxRowNormR1
% 
% vMinRowNormR1
%
% limits_t

lowerLimit = limits_t(1);
upperLimit = limits_t(2);

vec_x = lowerLimit:1:upperLimit;

vec_metric = vMinRowNormR1 ./ vMaxRowNormR1;


figure_name = '';
figure('name',figure_name)
hold on
plot(vec_x, log10(vec_metric))
hold off

end
