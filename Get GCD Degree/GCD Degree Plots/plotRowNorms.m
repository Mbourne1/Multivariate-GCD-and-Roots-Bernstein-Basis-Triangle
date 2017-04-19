function [] = plotRowNorms(arr_RowNorms, limits_k, limits_t)
%
% % Inputs
%
% arr_RowNorms : (Array)
%
% limits_k : (Int Int)
%
% limits_t : (Int Int)


lowerLimit_k = limits_k(1);
upperLimit_k = limits_k(2);

lowerLimit_t = limits_t(1);
upperLimit_t = limits_t(2);


figure_name = 'Plotting Row Norms';
figure('name',figure_name)
hold on

for i = lowerLimit_k:1:upperLimit_k
    
    vec = arr_RowNorms{i};
    
    vec_i = i.*ones(length(vec),1);
   
    plot(vec_i, vec); 
    
end

vline(lowerLimit_t);
vline(upperLimit_t);

hold off

end