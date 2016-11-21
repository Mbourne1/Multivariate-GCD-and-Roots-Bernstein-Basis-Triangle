function [] = o_gcd_batch()


% arr_ex_num = {'1','2','3','4','5','6','7','8','9','10'};
arr_ex_num = {'1','2','3','4'};
arr_el = {1e-8,1e-10};
arr_mean_method = {'Geometric Mean Matlab Method','None'};
arr_bool_alpha_theta = {'y','n'};
%arr_low_rank_approx_method = {'Standard SNTLN','Standard STLN','None'};
arr_low_rank_approx_method = {'None'};
arr_apf_method = {'None'};
arr_sylvester_matrix_type = {'T','DT','DTQ','TQ'};
%arr_sylvester_matrix_type = {'DTQ'};


parfor i1 = 1:1:length(arr_ex_num)
    for i2 = 1:1:length(arr_el)
        for i3 = 1:1:length(arr_mean_method)
            for i4 = 1:1:length(arr_bool_alpha_theta)
                for i5 = 1:1:length(arr_low_rank_approx_method)
                    for i6 = 1:1:length(arr_apf_method)
                        
                        for i7 = 1:1:length(arr_sylvester_matrix_type)
                            
                            ex_num = arr_ex_num{i1};
                            el = arr_el{i2};
                            em = 1e-12;
                            mean_method = arr_mean_method{i3};
                            bool_alpha_theta = arr_bool_alpha_theta{i4};
                            low_rank_approx_method = arr_low_rank_approx_method{i5};
                            apf_method = arr_apf_method{i6};
                            sylvester_matrix_type = arr_sylvester_matrix_type{i7};
                            
                            try
                                close all;
                                clc;
                                o_gcd(ex_num, el, em, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, sylvester_matrix_type)
                                fileId = fopen('log.txt','a')
                                fprintf(fileId,'%s','success \n');
                                fclose(fileId);
                            catch err
                                fileId = fopen('log.txt','a')
                                fprintf(fileId,'%s \n\n\n',getReport(err));
                                fclose(fileId);
                            end
                            
                        end
                    end
                end
            end
        end
    end
end


end