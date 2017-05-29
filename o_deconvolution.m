function [] = o_deconvolution(ex_num, emin, bool_preproc)
% Perform tests on deconvolution methods
%
%
% % Inputs
%
% ex_num : (String) Example number (String)
%
% bool_preproc : (Boolean) whether preprocessing is included.
%
% emin : (Float) Noise level
%
%
% % Outputs.
%
%
%
% % Example
%
% >> o_deconvolution('1', true, 1e-10);

% add location for Bivariate_Deconvolution_Examples() function
addpath(genpath('../Examples'));

% Add that folder plus all subfolders to the path.
addpath(genpath(pwd));

% Add symbolic variables
syms x y;

global SETTINGS
SETTINGS.SEED = 1024;
SETTINGS.PREPROC_DECONVOLUTIONS = bool_preproc;
SETTINGS.MAX_ERROR_DECONVOLUTIONS = 1e-14;
SETTINGS.MAX_ITERATIONS_DECONVOLUTIONS = 10;
SETTINGS.PLOT_GRAPHS = true;

% Get a matrix containing the symbolic factors of a polynomial and a vector
% containing the multiplicity of each factor.
[factor_mult_arr] = Deconvolution_Examples_Bivariate(ex_num);

factor = factor_mult_arr(:,1);

vMult = factor_mult_arr(:,2);
% vMult is symbolic so make it numeric
vMult = double(vMult);

% Get the highest power of any factor
highest_pwr = max(vMult);


% %
% %
% Generate the polynomials f_{0},...,f_{m}(x)

arr_sym_fxy = cell(highest_pwr + 1,1);

vDeg_arr_fxy = zeros(highest_pwr + 1,1);

for i = 0:1:highest_pwr
    
    % Get the multiplicities of the roots of f_{i+1}
    mults = ((vMult - i) + abs(vMult-i)) ./ 2;
    
    % Get the symbolic polynomial f_{i+1}
    arr_sym_fxy{i+1} = prod(factor.^(mults));
    
    % Get the degree of polynomial f_{i+1}(x)
    vDeg_arr_fxy(i+1) = double(feval(symengine, 'degree', (arr_sym_fxy{i+1})));
end

display(arr_sym_fxy{1})

% Get the degree structure of the polynomials h_{i} where h_{i} =
% f_{i-1}/f_{i}.
vDegt_arr_hxy = -1 .* diff(vDeg_arr_fxy);

% Get the sequence of polynomials h_{i}(x) in symbolic form
arr_sym_h = cell(length(arr_sym_fxy)-1, 1);
for i = 1:1:length(arr_sym_fxy) - 1
    
    arr_sym_h{i} = arr_sym_fxy{i} / arr_sym_fxy{i+1};
    
end


% %
% %
% Get coefficients vectors of f_{i}(x) and h_{i}(x)
nPolys_arr_fxy = length(arr_sym_fxy);
nPolys_arr_hxy = length(arr_sym_fxy) - 1;

arr_fxy_pwr = cell(nPolys_arr_fxy, 1);
arr_fxy_brn = cell(nPolys_arr_fxy, 1);

arr_hxy_pwr_exact = cell(nPolys_arr_hxy, 1);
arr_hxy_brn_exact = cell(nPolys_arr_hxy, 1);

for i = 1:1:nPolys_arr_fxy
    
    
    arr_fxy_pwr{i,1} = double(rot90(coeffs(arr_sym_fxy{i},[x,y],'All'),2));
    
    arr_fxy_brn{i,1} = PowerToBernstein(arr_fxy_pwr{i,1},vDeg_arr_fxy(i));
    
    if i <= nPolys_arr_hxy
        
        arr_hxy_pwr_exact{i,1} = double(rot90(coeffs(arr_sym_h{i},[x,y],'All'),2));
        
        arr_hxy_brn_exact{i,1} = PowerToBernstein(arr_hxy_pwr_exact{i,1},vDegt_arr_hxy(i));
        
    else
        arr_fxy_pwr{i,1} = 1;
    end
    
end

% %
% %
% %
% Add noise to the coefficients of f_{i}(x)
arr_fxy_brn_noisy = cell(nPolys_arr_fxy,1);

for i = 1:1:nPolys_arr_fxy
    
    arr_fxy_brn_noisy{i,1} = AddNoiseToPoly(arr_fxy_brn{i},emin);
    
end

% Define an array of deconvolution methods to be used
arr_DeconvolutionMethod = {...
    'Deconvolution Separate',...
    'Deconvolution Batch',...
    'Deconvolution Batch With STLN',...
    'Deconvolution Batch Constrained',...
    'Deconvolution Batch Constrained With STLN'};


nMethods = length(arr_DeconvolutionMethod);

% Testing deconvolution
LineBreakLarge();
arr_hx = cell(nMethods,1);
arr_Error = cell(nMethods,1);

for i = 1 : 1 : nMethods
    
    % Get deconvolution method
    method_name = arr_DeconvolutionMethod{i};
    
    switch method_name
        
        case 'Deconvolution Separate'
            arr_hx{i,1} = Deconvolve_Bivariate_Separate(arr_fxy_brn_noisy);
            
        case 'Deconvolution Batch'
            arr_hx{i,1} = Deconvolve_Bivariate_Batch(arr_fxy_brn_noisy, vDeg_arr_fxy);
            
        case 'Deconvolution Batch With STLN'
            arr_hx{i,1} = Deconvolve_Bivariate_Batch_With_STLN(arr_fxy_brn_noisy, vDeg_arr_fxy);
            
        case 'Deconvolution Batch Constrained'
            arr_hx{i,1} = Deconvolve_Bivariate_Batch_Constrained(arr_fxy_brn_noisy, vDeg_arr_fxy);
            
        case 'Deconvolution Batch Constrained With STLN'
            arr_hx{i,1} = Deconvolve_Bivariate_Batch_Constrained_With_STLN(arr_fxy_brn_noisy, vDeg_arr_fxy);
            
        otherwise
            error('err')
            
    end
    
    fprintf([mfilename ' : ' sprintf('%s \n', method_name )]);
    
    arr_Error{i,1} = GetErrors(arr_hx{i}, arr_hxy_brn_exact);
    
end





% Plotting


nPolys_hx = size(arr_hx,1);

if (SETTINGS.PLOT_GRAPHS)
    
    figure_name = sprintf([mfilename ':' 'Deconvolution Methods Error']);
    figure('name',figure_name)
    hold on
    
    for i = 1 : 1 : nMethods
        
        methodName = arr_DeconvolutionMethod{i};
        vError = arr_Error{i};
        plot(log10(vError), '-o', 'DisplayName', methodName)
        
    end
    
    
    
    legend(gca,'show');
    xlim([1 nPolys_hx]);
    xlabel('Factor')
    ylabel('log_{10} error')
    hold off
    
end


%--------------------------------------------------------------------------
% Console writing

for i = 1:1:nMethods
    
    methodName = arr_DeconvolutionMethod{i};
    vError = arr_Error{i};
    display([mfilename ' : ' sprintf('Error %s : %e', methodName, norm(vError))]);
    
end



% Initialise array to store error for each method
arr_ErrorNorm = cell(nMethods,1);

for i = 1 : 1 : nMethods
    arr_ErrorNorm{i} = norm(arr_Error{i});
end


PrintToResultsFile(ex_num, bool_preproc, emin, arr_DeconvolutionMethod, arr_ErrorNorm);


end

function [] = PrintToResultsFile(ex_num, bool_preproc, noise, arr_DeconvolutionMethod, arr_ErrorNorm)
% Print results of each deconvolution method to 
% 
% % Inputs
%
% ex_num : (String)
%
% bool_preproc : (Boolean)
%
% noise : (Float) 
%
% arr_DeconvolutionMethod : (Array of Strings)
%
% arr_ErrorNorm : (Array of Integers)
%



% Get number of deconvolution methods considered
nMethods = length(arr_DeconvolutionMethod);

fullFileName = sprintf('Deconvolution/Results/Results_o_deconvolutions%s.txt',datetime('today'));

% If file already exists append a line
if exist(fullFileName, 'file')
    
    fileID = fopen(fullFileName,'a');
    
    for i = 1 : 1 : nMethods
        
        method_name = arr_DeconvolutionMethod{i};
        error_norm = arr_ErrorNorm{i};
        WriteNewLine(method_name, error_norm)
        
    end
    fclose(fileID);
    
else % File doesnt exist so create it
    
    fileID = fopen( fullFileName, 'wt' );
    WriteHeader()
    
    for i = 1 : 1 : nMethods
        
        method_name = arr_DeconvolutionMethod{i};
        error_norm = arr_ErrorNorm{i};
        WriteNewLine(method_name, error_norm)
        
    end
    
    fclose(fileID);
    
end

    function WriteNewLine(method_name, error_norm)
        %
        % % Inputs
        %
        % method_name : (String)
        %
        % error_norm : (Float)
        %
        fprintf(fileID,'%s,%s,%s,%s,%s,%s \n',...
            datetime('now'),...
            ex_num,...
            num2str(bool_preproc),...
            num2str(noise),...
            method_name,...
            num2str(error_norm)...
        );
        
    end

    function WriteHeader()
        
        fprintf(fileID,'DATE, EX_NUM, BOOL_PREPROC, NOISE, Method, Error \n');
        
    end


end


function [vError] = GetErrors(arr_hxy_comp, arr_hxy_exact)

vError = zeros(length(arr_hxy_comp),1);

% Compare each computed h{i} with actual h_{i}
for i = 1:1:length(arr_hxy_comp)
    
    exact = arr_hxy_exact{i}./arr_hxy_exact{i}(1,1);
    comp = arr_hxy_comp{i}./ arr_hxy_comp{i}(1,1);
    
    vError(i) = norm(exact - comp) ./ norm(exact);
    
end

%display(vError);

end


function GetErr(vErrors,name)

err = norm(vErrors);
fprintf([mfilename ' : ' sprintf('Error in %s Method : %e',name,err) '\n'])

end