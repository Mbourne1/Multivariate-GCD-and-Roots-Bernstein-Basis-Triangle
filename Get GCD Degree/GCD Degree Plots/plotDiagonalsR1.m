function [] = plotDiagonalsR1(arr_R1, my_limits, limits_t)
%
% % Inputs
%
% arr_R1 : 
%
% my_limits :
%
% limits_t : [lower upper] :
%
% 

my_lowerLimit = my_limits(1);
my_upperLimit = my_limits(2);

lowerLimit = limits_t(1);
upperLimit = limits_t(2);

nSubresultants = my_upperLimit - my_lowerLimit + 1;

global SETTINGS

figure_name = sprintf('Diagonals of R1 from QR decomposition of %s',SETTINGS.SYLVESTER_BUILD_METHOD);
figure('name',figure_name)
hold on
for i = 1:1:nSubresultants
   
    % Get diagonal entries of R1
    vec = diag(arr_R1{i});
    
    % Get vector of [ i i i i ...]
    v_i = i*ones(length(vec));
    
    % plot
    plot(v_i,log10(vec),'*')
    
    
end

% Add vertical lines to show limits
vline(lowerLimit);
vline(upperLimit);
hold off

end