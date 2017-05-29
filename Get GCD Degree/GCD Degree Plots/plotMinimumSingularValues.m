function [] = plotMinimumSingularValues(vMinimumSingularValues, limits_k, limits_t, rank_range)
%
% % Inputs
%
% vMinimumSingularValues : (Vector)
%
% limits_k : (Int Int)
%
% limits_t : (Int Int)
%
% rank_range : [Float Float]

%
lowerLimit_k = limits_k(1);
upperLimit_k = limits_k(2);

%
lowerLimit_t = limits_t(1);
upperLimit_t = limits_t(2);

%
rank_range_low = rank_range(1);
rank_range_high = rank_range(2);

vec_x = lowerLimit_k : 1 : upperLimit_k;

figure_name = sprintf('');
figure('name',figure_name)
hold on
plot(vec_x, log10(vMinimumSingularValues),'-s')


% 
hline(rank_range_low);
hline(rank_range_high);
hline(mean(rank_range));

%
vline(lowerLimit_t);
vline(upperLimit_t);


hold off

end