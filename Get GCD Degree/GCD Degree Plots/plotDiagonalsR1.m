function [] = plotDiagonalsR1(arr_R1, myLimits, limits_t)
%
% % Inputs
%
% arr_R1 : 
%
% myLimits :
%
% limits_t : [lower upper] :
%
% 

% Get my limits
myLowerLimit = myLimits(1);
myUpperLimit = myLimits(2);

% Get actual limits
lowerLimit = limits_t(1);
upperLimit = limits_t(2);

% Get number of subresultants constructed
nSubresultants = myUpperLimit - myLowerLimit + 1;

global SETTINGS

figure_name = sprintf('Diagonals of R1 from QR decomposition of %s',SETTINGS.SYLVESTER_BUILD_METHOD);
figure('name',figure_name)
hold on

for i = 1 : 1 : nSubresultants
   
    k = myLowerLimit + (i-1);
    
    % Get diagonal entries of R1
    temp_vec = diag(arr_R1{i});
    
    % Get vector of [ i i i i ...]
    vec_x = k*ones(length(temp_vec));
    
    % plot
    plot(vec_x, log10(temp_vec), '*')
    
    
end

% Add vertical lines to show limits
vline(lowerLimit);
vline(upperLimit);
hold off

end