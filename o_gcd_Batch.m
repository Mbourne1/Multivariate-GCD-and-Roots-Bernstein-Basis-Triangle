function [] = o_Batch


arr_ex_num = {'1','2'};
arr_el = {1e-10};
arr_mean_method = {'Geometric Mean Matlab Method','None'};
arr_bool_alpha_theta = {'y','n'};
arr_low_rank_approx_method = {'Standard STLN','None'};

for i1 = 1:1:length(arr_ex_num)
    for i2 = 1:1:length(arr_el)
        for i3 = 1:1:length(arr_mean_method)
            for i4 = 1:1:length(arr_bool_alpha_theta)
                for i5 = 1:1:length(arr_low_rank_approx_method)
                    
                    ex_num = arr_ex_num{i1};
                    el = arr_el{i2};
                    mean_method = arr_mean_method{i3};
                    bool_alpha_theta = arr_bool_alpha_theta{i4};
                    low_rank_approx_method = arr_low_rank_approx_method{i5};
                    
                    o_gcd(ex_num,el,mean_method,bool_alpha_theta,low_rank_approx_method)
                end
            end
        end
    end
end


end