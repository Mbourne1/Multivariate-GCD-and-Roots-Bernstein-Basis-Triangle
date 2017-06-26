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

vec_x = lowerLimit_k : 1 : upperLimit_k;

figure_name = sprintf('');
figure('name',figure_name)
hold on
plot(vec_x, log10(vMinimumSingularValues),'-s')

xlim(limits_k);
xlabel('k')
ylabel('log_{10}(\sigma_{k})')
% 
hline(rank_range,{'r','r'});
vline(limits_t, {'r','r'});



hold off

end