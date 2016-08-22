function [] = Test_Deconvolution(ex_num,noise,bool_preproc)
% Perform tests on deconvolution methods
%
%
% % Inputs
%
% ex_num : Example number (String)
%
% emin : Noise level
%
%
% % Outputs.
%
%  
%
% Example
%
% >> Test_Deconvolution('1',1e-10);



syms x y;

global SETTINGS
SETTINGS.SEED = 1024;
SETTINGS.PREPROC_DECONVOLUTIONS = bool_preproc;
SETTINGS.MAX_ERROR_DECONVOLUTIONS = 1e-14;
SETTINGS.MAX_ITERATIONS_DECONVOLUTIONS = 10;

switch ex_num
    case '1' 
        
        % (x + y + 0.5)^7 * (x+y+1.2)^12
        
        factor(1,1) = (x + y + 0.5);
        factor(2,1) = (x + y + 1.2);
        
        % Set the multiplicity of each factor
        vMult = [7; 12];
        

    case '2'
        factor(1,1) = (x + y + 0.5);
        factor(2,1) = (x + y + 1.2);
        
        % Set the multiplicity of each factor
        vMult = [10; 15];
    
    case '3'
        % (x + y + 0.5)^7(x + y + 1.2)^10(x + y - 0.15)^15
        factor(1,1) = (x + y + 0.5);
        factor(2,1) = (x + y + 1.2);
        factor(3,1) = (x + y - 0.15);
        
        
        % Set the multiplicity of each factor
        vMult = [7; 10; 15];
        
end


% Get the highest power of any factor
highest_pwr = max(vMult);


% %
% %
% Generate the polynomials f_{0},...,f_{m}(x)

arr_sym_fxy = cell(highest_pwr+1,1);
vDegt_arr_fxy = zeros(highest_pwr+1,1);
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
vDegt_arr_hxy = -1.*diff(vDegt_arr_fxy);

% Get the degree structure of the polynomials w_{i}
vDegt_arr_wxy = diff([vDegt_arr_hxy; 0]);

% Get the multiplicities of the roots
vMultiplicities = find(vDegt_arr_wxy~=0);

% Get the sequence of polynomials h_{i}(x) in symbolic form
arr_sym_h = cell(length(arr_sym_fxy)-1,1);
for i = 1:1:length(arr_sym_fxy)-1
    arr_sym_h{i} = arr_sym_fxy{i} / arr_sym_fxy{i+1};
end


% %
% %
% Get coefficients vectors of f_{i}(x) and h_{i}(x)
nPolys_arr_fxy = length(arr_sym_fxy);
nPolys_arr_hxy = length(arr_sym_fxy) - 1;

arr_fxy_pwr = cell(nPolys_arr_fxy,1);
arr_fxy_brn = cell(nPolys_arr_fxy,1);

arr_hxy_pwr = cell(nPolys_arr_hxy,1);
arr_hxy_brn = cell(nPolys_arr_hxy,1);

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
    arr_fxy_brn_noisy{i,1} = Noise(arr_fxy_brn{i},noise);
end

%--------------------------------------------------------------------------
% %
% %
% %
% Testing deconvolution
LineBreakLarge();
fprintf([mfilename ' : ' 'Separate \n'])

arr_hxy_brn_separate = Deconvolve_Bivariate_Separate(arr_fxy_brn_noisy);
vErrors_separate = GetErrors(arr_hxy_brn_separate,arr_hxy_brn);

% -------------------------------------------------------------------------
% %
% %
% Testing deconvolution batch method with constraints and low rank approx
LineBreakLarge()
fprintf([mfilename ' : ' 'Batch \n'])

arr_hxy_brn_batch = Deconvolve_Bivariate_Batch(arr_fxy_brn_noisy,vDegt_arr_fxy);
vErrors_batch = GetErrors(arr_hxy_brn_batch,arr_hxy_brn);

% -------------------------------------------------------------------------
% %
% %
% Testing deconvolution batch method with constraints and low rank approx
LineBreakLarge()
fprintf([mfilename ' : ' 'Batch with STLN \n'])

arr_hxy_brn_batch_with_STLN = Deconvolve_Bivariate_Batch_With_STLN(arr_fxy_brn_noisy,vDegt_arr_fxy);
vErrors_batch_with_STLN = GetErrors(arr_hxy_brn_batch_with_STLN,arr_hxy_brn);




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


fileID = fopen('Deconvolution/Test_Deconvolution.txt','a');
fprintf(fileID,'%s %s %6.2e %6.2e %6.2e %6.2e %6.2e %6.2e\n',ex_num,bool_preproc,noise,A);
fclose(fileID);

end


function [vError] = GetErrors(arr_hxy_comp, arr_hxy_exact)

vError = zeros(length(arr_hxy_comp),1);

% Compare each computed h{i} with actual h_{i}
for i = 1:1:length(arr_hxy_comp)
    
    exact = arr_hxy_exact{i}./arr_hxy_exact{i}(1,1);
    comp = arr_hxy_comp{i}./ arr_hxy_comp{i}(1,1);
    
    vError(i) = norm(exact - comp) ./ norm(exact);
    
end

display(vError);

end


function GetErr(vErrors,name)

err = norm(vErrors);
fprintf([mfilename ' : ' sprintf('Error in %s Method : %e',name,err) '\n'])

end