function [] = plotRowNorms(arr_RowNorms, limits_t)
%
% % Inputs
%
% arr_RowNorms : 
%
% limits_t :


lowerLimit = limits_t(1);
upperLimit = limits_t(2);


figure_name = 'Plotting Row Norms';
figure('name',figure_name)
hold on

for i = lowerLimit:1:upperLimit
    
    vec = arr_RowNorms{i};
    
    vec_i = i.*ones(length(vec),1);
   
    plot(vec_i, vec); 
    
end

hold off

end