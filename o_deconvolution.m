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

% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 

% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

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

vDegt_arr_fxy = zeros(highest_pwr + 1,1);

for i = 0:1:highest_pwr
    
    % Get the multiplicities of the roots of f_{i+1}
    mults = ((vMult - i) + abs(vMult-i)) ./ 2;
    
    % Get the symbolic polynomial f_{i+1}
    arr_sym_fxy{i+1} = prod(factor.^(mults));
    
    % Get the degree of polynomial f_{i+1}(x)
    vDegt_arr_fxy(i+1) = double(feval(symengine, 'degree', (arr_sym_fxy{i+1})));
end

display(arr_sym_fxy{1})

% Get the degree structure of the polynomials h_{i} where h_{i} =
% f_{i-1}/f_{i}.
vDegt_arr_hxy = -1 .* diff(vDegt_arr_fxy);

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

arr_hxy_pwr = cell(nPolys_arr_hxy, 1);
arr_hxy_brn = cell(nPolys_arr_hxy, 1);

for i = 1:1:nPolys_arr_fxy
    
    
    arr_fxy_pwr{i,1} = double(rot90(coeffs(arr_sym_fxy{i},[x,y],'All'),2));
    
    arr_fxy_brn{i,1} = PowerToBernstein(arr_fxy_pwr{i,1},vDegt_arr_fxy(i));
    
    if i <= nPolys_arr_hxy
        
        arr_hxy_pwr{i,1} = double(rot90(coeffs(arr_sym_h{i},[x,y],'All'),2));
        
        arr_hxy_brn{i,1} = PowerToBernstein(arr_hxy_pwr{i,1},vDegt_arr_hxy(i));
        
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

%--------------------------------------------------------------------------
% %
% %
% %
% Testing deconvolution separate method
LineBreakLarge();
fprintf([mfilename ' : ' 'Separate \n'])

arr_hxy_brn_separate = Deconvolve_Bivariate_Separate(arr_fxy_brn_noisy);
vErrors_separate = GetErrors(arr_hxy_brn_separate, arr_hxy_brn);

% -------------------------------------------------------------------------
% %
% %
% Testing deconvolution batch method with constraints and low rank approx
LineBreakLarge()
fprintf([mfilename ' : ' 'Batch \n'])

arr_hxy_brn_batch = Deconvolve_Bivariate_Batch(arr_fxy_brn_noisy, vDegt_arr_fxy);
vErrors_batch = GetErrors(arr_hxy_brn_batch, arr_hxy_brn);

% -------------------------------------------------------------------------
% %
% %
% Testing deconvolution batch method with constraints and low rank approx
LineBreakLarge()
fprintf([mfilename ' : ' 'Batch with STLN \n'])

arr_hxy_brn_batch_with_STLN = Deconvolve_Bivariate_Batch_With_STLN(arr_fxy_brn_noisy,vDegt_arr_fxy);
vErrors_batch_with_STLN = GetErrors(arr_hxy_brn_batch_with_STLN, arr_hxy_brn);

% -------------------------------------------------------------------------
% %
% %
% %
% Testing deconvolution batch method with constraints
LineBreakLarge()
fprintf([mfilename ' : ' 'Batch Constrained \n'])

arr_hxy_brn_batch_constrained = Deconvolve_Bivariate_Batch_Constrained(arr_fxy_brn_noisy,vDegt_arr_fxy);
vErrors_batch_constrained = GetErrors(arr_hxy_brn_batch_constrained,arr_hxy_brn);

% -------------------------------------------------------------------------
% %
% %
% %
% Testing deconvolution batch method with constraints
LineBreakLarge()
fprintf([mfilename ' : ' 'Batch Constrained with STLN \n'])

arr_hxy_brn_batch_constrained_with_STLN = Deconvolve_Bivariate_Batch_Constrained_With_STLN(arr_fxy_brn_noisy,vDegt_arr_fxy);
vErrors_batch_constrained_with_STLN = GetErrors(arr_hxy_brn_batch_constrained_with_STLN,arr_hxy_brn);




%-------------------------------------------------------------------------
% %
% %
% %
% %
% Plotting

if( SETTINGS.PLOT_GRAPHS)
    
    
    
    fig_name = sprintf([mfilename ' : ' 'Plotting Errors']);
    figure('name',fig_name)
    hold on
    plot(log10(vErrors_separate),'-s','DisplayName','Separate')
    plot(log10(vErrors_batch),'-*','DisplayName','Batch')
    plot(log10(vErrors_batch_with_STLN),'-*','DisplayName','Batch STLN')
    plot(log10(vErrors_batch_constrained),'-o','DisplayName','Batch Constrained')
    plot(log10(vErrors_batch_constrained_with_STLN),'-o','DisplayName','Batch Constrained STLN')
    legend(gca,'show')
    xlabel('i')
    ylabel('log_{10} Error in h_{i}(x,y)')
    hold off
    
end
%--------------------------------------------------------------------------
GetErr(vErrors_separate,'Separate')
GetErr(vErrors_batch,'Batch')
GetErr(vErrors_batch_with_STLN,'Batch With STLN')
GetErr(vErrors_batch_constrained,'Batch Constrained')
GetErr(vErrors_batch_constrained_with_STLN,'Batch Constrained STLN');
%--------------------------------------------------------------------------

A = ...
    [
    norm(vErrors_separate),...
    norm(vErrors_batch),...
    norm(vErrors_batch_with_STLN),...
    norm(vErrors_batch_constrained),...
    norm(vErrors_batch_constrained_with_STLN)
    ];

my_error.separate = norm(vErrors_separate);
my_error.batch = norm(vErrors_batch);
my_error.batchSTLN = norm(vErrors_batch_with_STLN);
my_error.batchConstrained = norm(vErrors_batch_constrained);
my_error.batchConstrainedSTLN = norm(vErrors_batch_constrained_with_STLN);

PrintToResultsFile(ex_num,bool_preproc,emin,my_error);

end

function [] = PrintToResultsFile(ex_num,bool_preproc,noise,my_error)


fullFileName = sprintf('Deconvolution/Results/Results_o_deconvolutions%s.txt',datetime('today'));

% If file already exists append a line
if exist(fullFileName, 'file')
    
    fileID = fopen(fullFileName,'a');
    WriteNewLine()
    fclose(fileID);
    
else % File doesnt exist so create it
    
    fileID = fopen( fullFileName, 'wt' );
    WriteHeader()
    WriteNewLine()
    fclose(fileID);
    
end

    function WriteNewLine()
        
        %
        fprintf(fileID,'%s,%s,%s,%s,%s,%s,%s,%s,%s \n',...
            datetime('now'),...
            ex_num,...
            bool_preproc,...
            noise,...
            my_error.separate,...
            my_error.batch,...
            my_error.batchSTLN ,...
            my_error.batchConstrained ,...
            my_error.batchConstrainedSTLN....
            );
        
    end

    function WriteHeader()
        
        fprintf(fileID,'DATE,EX_NUM,BOOL_PREPROC,NOISE,err_separate,err_batch,err_batch_STLN,err_constrained,err_constrained_STLN \n');
        
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