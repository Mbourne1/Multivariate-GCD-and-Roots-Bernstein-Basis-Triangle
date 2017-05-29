function [] = plotMaxMinDiagonalR1(vRatio_MaxMin_DiagonalEntry, limits_k, limits_t, rank_range)
%
% % Inputs
%
% vRatio_MaxMin_DiagonalEntry : (Vector)
%
% limits_k : (Int Int)
%
% limits_t : (Int Int)
%
% rank_range : [Float Float]

global SETTINGS

%
lowerLimit_k = limits_k(1);
upperLimit_k = limits_k(2);

%
lowerLimit_t = limits_t(1);
upperLimit_t = limits_t(2);

%
lower_rank_range = rank_range(1);
upper_rank_range = rank_range(2);

figure_name = [mfilename ' : ' sprintf('Max;min Diagonals of R1 from QR Decomposition of %s', SETTINGS.SYLVESTER_BUILD_METHOD)];
    
figure('name', figure_name)
hold on
x_vec = lowerLimit_k : 1 : upperLimit_k;
plot(x_vec, log10(vRatio_MaxMin_DiagonalEntry),'-s')

%
hline(lower_rank_range,'-r');
hline(upper_rank_range,'-r');
hline(mean(rank_range),'-b');

%
vline(lowerLimit_t);
vline(upperLimit_t);


hold off


end