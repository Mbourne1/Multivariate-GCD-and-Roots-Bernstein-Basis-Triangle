function [] = plotMinimumSingularValues(vMinimumSingularValues,myLimits, limits_t)
%
% % Inputs
%
% vMinimumSingularValues :
%
% myLimits :
%
% limits_t :


lowerLimit = myLimits(1);
upperLimit = myLimits(2);

lowerLimit_t = limits_t(1);
upperLimit_t = limits_t(2);

vec_x = lowerLimit : 1 : upperLimit;

figure_name = sprintf('');
figure('name',figure_name)
hold on
plot(vec_x, log10(vMinimumSingularValues),'-s')

vline(lowerLimit_t);
vline(upperLimit_t);


hold off

end