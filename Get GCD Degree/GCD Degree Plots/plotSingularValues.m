 function [] = plotSingularValues(arr_SingularValues, limits_k, limits_t)
%
% % Inputs
%
% arr_SingularValues : (Array)
%
% myLimits :
%
% limits_t :

lowerLimit_k = limits_k(1);
upperLimit_k = limits_k(2);

lowerLimit_t = limits_t(1);
upperLimit_t = limits_t(2);

nSubresultants = upperLimit_k - lowerLimit_k + 1;

figure_name = 'Singular Values';
figure('name',figure_name)



hold on
for i = 1:1:nSubresultants
   
    k = lowerLimit_k + (i-1);
    
    temp_vec = arr_SingularValues{i};
    
    vec_k = k.*ones(length(temp_vec));
    
    plot(vec_k,log10(temp_vec),'*');
    
end

vline(lowerLimit_t);
vline(upperLimit_t);

hold off

end