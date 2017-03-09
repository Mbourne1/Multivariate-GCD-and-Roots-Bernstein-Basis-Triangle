function [] = o_gcd_Bivariate_2Polys_Batch
% Peform multiple gcd computations


% Initialise arrays
arr_ex_num = {'1', '2','3','4'};

arr_el = {1e-12,1e-10,1e-8,1e-6};

arr_mean_method = ...
    {...
    'Geometric Mean Matlab Method',...
    'Geometric Mean My Method',...
    'None'...
    };

arr_bool_alpha_theta = {true,false};

arr_low_rank_approx_method = {'None','Standard STLN'};

parpool(length(2));



parfor i1 = 1:1:length(arr_ex_num)
    for i2 = 1:1:length(arr_el)
        for i3 = 1:1:length(arr_mean_method)
            for i4 = 1:1:length(arr_bool_alpha_theta)
                for i5 = 1:1:length(arr_low_rank_approx_method)

                    ex_num = arr_ex_num{i1};
                    el = arr_el{i2};
                    mean_method = arr_mean_method{i3};
                    bool_alpha_theta = arr_bool_alpha_theta{i4};
                    low_rank_approx_method = arr_low_rank_approx_method{i5};
                    
         
                    try

                        o_gcd_Bivariate_2Polys(ex_num, el, mean_method, bool_alpha_theta, low_rank_approx_method);

                    catch error
                        error.message
                        WriteToErrorFile(ex_num,el,mean_method,bool_alpha_theta)
                    end
                end
            end
        end
    end
end




end

function []= WriteToErrorFile(ex_num,noise,mean_method,bool_alpha_theta)



fullFileName = 'Results_GCD_errors.txt';


if exist(fullFileName, 'file')
    fileID = fopen(fullFileName,'a');
    fprintf(fileID,'%s, \t %s, \t %s, \t %s, \t %s \n',...
        datetime('now'),...
        ex_num,...
        num2str(noise),...
        mean_method,...
        bool_alpha_theta...
        );
    fclose(fileID);
else
    % File does not exist.
    warningMessage = sprintf('Warning: file does not exist:\n%s', fullFileName);
    uiwait(msgbox(warningMessage));
end
end