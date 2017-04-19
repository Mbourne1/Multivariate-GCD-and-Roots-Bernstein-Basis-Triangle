function [] = o_gcd_Bivariate_2Polys_Batch
% Peform multiple gcd computations


% Initialise arrays
arr_ex_num = {'1', '2','3','4'};

arr_el = {1e-12,1e-10,1e-8,1e-6};

arr_mean_method = ...
    {...
    'Geometric Mean Matlab Method',...
    'None'...
    };

arr_bool_alpha_theta = {true,false};

arr_low_rank_approx_method = {'None','Standard STLN'};

%parpool(length(2));

apf_method = 'None';
sylvester_build_method = 'DTQ';


for i1 = 1:1:length(arr_ex_num)
    for i2 = 1:1:length(arr_el)
        for i3 = 1:1:length(arr_mean_method)
            for i4 = 1:1:length(arr_bool_alpha_theta)
                for i5 = 1:1:length(arr_low_rank_approx_method)

                    ex_num = arr_ex_num{i1};
                    emin = arr_el{i2};
                    emax = 1e-10;
                    mean_method = arr_mean_method{i3};
                    bool_alpha_theta = arr_bool_alpha_theta{i4};
                    low_rank_approx_method = arr_low_rank_approx_method{i5};
                    
         
                     
                    
                    try

                        
                        close all;
                        clc;
                        o_gcd_Bivariate_2Polys(ex_num, emin, emax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, sylvester_build_method);
                        myFileName = 'log_GCD_Bivariate_2Polys.txt';
                        fileId = fopen(myFileName,'a');
                        fprintf(fileId,'%s %s \n',datetime('now'), 'success \n');
                        fclose(fileId);

                    catch err
                    
                        myFileName = 'log_GCD_Bivariate_2Polys.txt';
                        fileId = fopen(myFileName,'a');
                        fprintf(fileId,'%s %s \n\n\n', datetime('now'), getReport(err));
                        fclose(fileId);
                        
                    end
                end
            end
        end
    end
end




end