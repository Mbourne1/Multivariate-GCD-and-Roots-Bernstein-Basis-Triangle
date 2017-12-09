function [] = Experiment1SylvesterFormat_2Polys(ex_num, bool_preproc)
%
%
%
% % Inputs
%
% ex_num : (String) Example number
%
% Experiment1SylvesterFormat_2Polys('1')

addpath(genpath(pwd));

close all;
clc;

% Examples
% 1 : Degree too low, not complex enough - even high level noise, method
% works well for all formats
%
% 2 :
%
%
% Bad Examples, 1 2 3 4 5 6 8 9

% Good Examples

% ex_num = '20'; el = 1e-5; eu = 1e-5;






% Constants ---------------------------------------------------------------

el = 1e-6;
eu = 1e-5;


switch bool_preproc
    case true
        mean_method = 'Geometric Mean Matlab Method';
        bool_alpha_theta = true;
    case false
        mean_method = 'None';
        bool_alpha_theta = false;
end


%mean_method = 'Geometric Mean Matlab Method';


low_rank_approx_method = 'None';
apf_method = 'None';
rank_revealing_metric = 'Minimum Singular Values';



arrSylvesterFormat = {'T', 'DT', 'TQ', 'DTQ', 'DTQ Denominator Removed'};


for i = 1 : 1 : length(arrSylvesterFormat)
    
    
    sylvester_format = arrSylvesterFormat{i};
    try
        %close all; clc;
        o_gcd_Bivariate_2Polys(ex_num, el, eu, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, sylvester_format, rank_revealing_metric)
        %SavePlots()
        
    catch
    end
end



    function [] = SavePlots()
        %
        % % Inputs
        %
        % ex_num : (String)
        %
        % str : (String)
        %
        % sylvester_build_method
        
        directory_name = strcat('UnivariateSylvesterFormatFigures/Example',(ex_num),'/Figures/');
        
        mkdir(directory_name)
        h = get(0,'children');
        
        myFileName = strcat(sylvester_format);
        
        for j = 2 : 1 : 2
            try
                %saveas(h(i), [directory_name num2str(length(h) + 1 - i)], 'fig');
                saveas(h(j), [directory_name myFileName], 'fig');
                saveas(h(j), [directory_name myFileName], 'eps');
                saveas(h(j), [directory_name myFileName], 'png');
            catch
                error('err')
            end
        end
        
        
    end


end
