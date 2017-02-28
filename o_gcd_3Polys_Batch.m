function [] = o_gcd_3Polys_Batch
%
%
% Examples
%
% >> o_gcd_3Polys_Batch

% Peform multiple gcd computations


% Initialise arrays
arr_ex_num = {'1', '2','3','4'};

arr_emax = {1e-12, 1e-10, 1e-8, 1e-6};
emin = 1e-14;

arr_mean_method = ...
    {...
    'Geometric Mean Matlab Method'
    };

arr_bool_alpha_theta = {'y','n'};

arr_low_rank_approx_method = {'None'};
apf_method = {'None'};


parpool(length(2));



parfor i1 = 1:1:length(arr_ex_num)
    
    for i2 = 1:1:length(arr_emax)
        
        
        for i3 = 1:1:length(arr_mean_method)
            for i4 = 1:1:length(arr_bool_alpha_theta)
                for i5 = 1:1:length(arr_low_rank_approx_method)

                    ex_num = arr_ex_num{i1};
                    emax = arr_emax{i2};
                    mean_method = arr_mean_method{i3};
                    bool_alpha_theta = arr_bool_alpha_theta{i4};
                    low_rank_approx_method = arr_low_rank_approx_method{i5};
                    
         
                    try

                        
                        o_gcd_3Polys(ex_num, emin, emax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, sylvester_type)
                        
                        
                    catch error
                        error.message
                        WriteToErrorFile(ex_num,emax,emin,mean_method,bool_alpha_theta)
                    end
                end
            end
        end
    end
end




end

function []= WriteToErrorFile(ex_num, emax, emin, mean_method, bool_alpha_theta)



fullFileName = 'Results_GCD_errors.txt';


if exist(fullFileName, 'file')
    fileID = fopen(fullFileName,'a');
    fprintf(fileID,'%s, \t %s, \t %s, \t %s, \t %s \n',...
        datetime('now'),...
        ex_num,...
        num2str(emax),...
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