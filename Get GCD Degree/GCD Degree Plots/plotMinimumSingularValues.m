function [] = plotMinimumSingularValues(vMinimumSingularValues, limits_k, limits_t)
%
% % Inputs
%
% vMinimumSingularValues :
%
% limits_k : (Int Int)
%
% limits_t : (Int Int)


lowerLimit_k = limits_k(1);
upperLimit_k = limits_k(2);

lowerLimit_t = limits_t(1);
upperLimit_t = limits_t(2);

vec_x = lowerLimit_k : 1 : upperLimit_k;

figure_name = sprintf('');
figure('name',figure_name)
hold on
plot(vec_x, log10(vMinimumSingularValues),'-s')

vline(lowerLimit_t);
vline(upperLimit_t);


hold off

end