function [] = plotMaxMinDiagonalR1(vRatio_MaxMin_DiagonalEntry, myLimits, limits_t)
%
% % Inputs
%
% vRatio_MaxMin_DiagonalEntry
%
% myLimits
%
% limits_t


global SETTINGS

% Get my working limits
lowerLimit = myLimits(1);
upperLimit = myLimits(2);

lower_limit_t = limits_t(1);
upper_limit_t = limits_t(2);

figure_name = [mfilename ' : ' sprintf('Max;min Diagonals of R1 from QR Decomposition of %s', SETTINGS.SYLVESTER_BUILD_METHOD)];
    
figure('name',figure_name)
hold on
x_vec = lowerLimit : 1 : upperLimit;
plot(x_vec,log10(vRatio_MaxMin_DiagonalEntry),'-s')
vline(lower_limit_t);
vline(upper_limit_t);

hold off


end