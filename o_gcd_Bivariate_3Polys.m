function [] = o_gcd_Bivariate_3Polys(ex_num, emin, emax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, sylvester_type)
% o_gcd(ex_num,emin,mean_method,bool_alpha_theta,low_rank_approx_method,apf_method,degree_method)
%
% Calculate the GCD d(x,y) of two polynomials f(x,y) and g(x,y) taken from
% the example file.
%
% % Inputs.
%
% ex_num  : (String) Example Number (String)
%
% emin : (Float) Minimum Noise level
%
% emax : (Float) Maximum signal to noise ratio
%
% mean_method : (String)
%       'Geometric Mean Matlab Method'
%       'None'
%
% bool_alpha_theta (Boolean)
%       true - Include Preprocessing
%       false - Exclude Preprocessing
%
% low_rank_approx_method (String)
%       'Standard SNTLN'
%       'None'
%
% apf_method (String)
%       'None'
%       'Standard APF Linear'
%       'Standard APF Nonlinear'
%
%
% sylvester_matrix_type : (String)
%       * T
%       * DT
%       * DTQ
%       * TQ
%
% % Example
%
% >> o_gcd_Bivariate_3Polys('1', 1e-10, 1e-12, 'Geometric Mean Matlab Method', true, 'None', 'None', 'DTQ')
% >> o_gcd_Bivariate_3Polys('1', 1e-10, 1e-12, 'Geometric Mean Matlab Method', true, 'Standard STLN', 'Standard Nonlinear APF', 'DTQ')


% Set the Global Variables
global SETTINGS

% add path
restoredefaultpath

addpath(...
    'APF',...
    'Bernstein Functions',...
    'Build Matrices',...
    'Build Sylvester Matrix',...
    'Formatting',...
    'GCD Methods',...
    'Get Cofactors',...
    'Get GCD Coefficients',...
    'Results'...
    );

addpath(genpath('Examples'));
addpath(genpath('Get GCD Degree'));
addpath(genpath('Preprocessing'));
addpath(genpath('Low Rank Approximation'));


problem_type = 'GCD';

% % Ensure that minimum noise level is less than maximum noise level
if emin > emax
    temp = emin;
    emin = emax;
    emax = temp;
end

% Set global variables
SetGlobalVariables(problem_type, ex_num, emin, emax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, sylvester_type)

% Print Parameters to screen
fprintf('INPUTS. \n')
fprintf('EXAMPLE NUMBER %s \n',ex_num)
fprintf('EMIN : %s \n',emin)
fprintf('EMAX : %s \n',emax)
fprintf('MEAN METHOD : %s \n', mean_method)
fprintf('PREPROCESSING : %s \n',bool_alpha_theta)
fprintf('LOW RANK METHOD : %s \n',low_rank_approx_method)
fprintf('APF METHOD : %s \n', apf_method)
fprintf('SYLVESTER TYPE : %s \n', sylvester_type)

% Get example polynomials
[fxy_exact, gxy_exact, hxy_exact,...
    uxy_exact,vxy_exact, wxy_exact,...
    dxy_exact,...
    m,~,~,n,~,~,o,~,~,...
    ~, ~, ~] = Examples_GCD_3Polys(ex_num);




DisplayDegreeStructure_3Polys();

% %
% %
% Add Noise

% Add noise to the coefficients of f and g
[fxy, ~] = AddVariableNoiseToPoly(fxy_exact, emin, emax);
[gxy, ~] = AddVariableNoiseToPoly(gxy_exact, emin, emax);
[hxy, ~] = AddVariableNoiseToPoly(hxy_exact, emin, emax);

% %
% % Get the GCD by zengs method
%[u,v,w] = o_gcd_zeng(fxy,gxy);


% Get GCD d(x,y) and quotient polynomials u(x,y) and v(x,y)
lowerLimit = 1;
upperLimit = min([m,n,o]);
t_limits = [lower_limit, upperLimit];



switch SETTINGS.DEGREE_METHOD
    case 'Total'
        
        fxy_matrix_padd = zeros(m+1,m+1);
        gxy_matrix_padd = zeros(n+1,n+1);
        
        [m1, m2] = GetDegree_Bivariate(fxy);
        fxy_matrix_padd(1:m1+1, 1:m2+1) = fxy;
        
        [n1, n2] = size(gxy);
        gxy_matrix_padd(1:n1+1, 1:n2+1) = gxy;
        
        fxy = fxy_matrix_padd;
        gxy = gxy_matrix_padd;
        
    case 'Relative'
        
    case 'Both'
end

% Get the GCD by my method
[fxy, gxy, hxy, dxy, uxy, vxy, wxy, t] = ...
    o_gcd_mymethod_3Polys(fxy, gxy, hxy, m, n, o, t_limits);



% Get Error in u(x,y), v(x,y) and d(x,y)
my_error.uxy = GetError(uxy, uxy_exact);
my_error.vxy = GetError(vxy, vxy_exact);
my_error.wxy = GetError(wxy, wxy_exact);

my_error.dxy = GetError(dxy, dxy_exact);

% Print the error in u(x,y), v(x,y) and d(x,y)
fprintf([mfilename ' : ' sprintf('Distance u(x,y) : %e \n', my_error.uxy)]);
fprintf([mfilename ' : ' sprintf('Distance v(x,y) : %e \n', my_error.vxy)]);
fprintf([mfilename ' : ' sprintf('Distance w(x,y) : %e \n', my_error.wxy)]);
fprintf([mfilename ' : ' sprintf('Distance d(x,y) : %e \n', my_error.dxy)]);

PrintToResultsFile(m, n, o, t, my_error)

end

function [dist] = GetDistance(fxy,gxy)
%
% % Inputs
%
% fxy : (Matrix)
%
% gxy : (Matrix) 


fxy = fxy./fxy(1,1);
gxy = gxy./gxy(1,1);

try
    dist = norm(fxy - gxy) ./ norm(fxy);
catch
    dist = 1000;
end

end

function []= PrintToResultsFile(m,n,o,t,my_error)
%
% % Inputs
%
% m : (Int)
%
% n : (Int) 
%
% o : (Int)
%
% t : (Int)
%
% my_error :

global SETTINGS

fullFileName = sprintf('Results/Results_o_gcd_3Polys%s.txt',datetime('today'));

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
        
        % 19 FIELDS
        fprintf(fileID,'%s,%s,%s, %s,%s,%s, %s,%s,%s, %s,%s,%s, %s,%s,%s, %s,%s,%s, %s \n',...
            datetime('now'),...
            SETTINGS.EX_NUM,...
            int2str(m),...
            int2str(n),...
            int2str(o),...
            int2str(t),...
            my_error.uxy,...
            my_error.vxy,...
            my_error.wxy,...
            my_error.dxy,...
            SETTINGS.MEAN_METHOD,...
            SETTINGS.BOOL_ALPHA_THETA,...
            SETTINGS.LOW_RANK_APPROX_METHOD,...
            num2str(SETTINGS.LOW_RANK_APPROX_REQ_ITE),...
            SETTINGS.APF_METHOD,...
            num2str(SETTINGS.APF_REQ_ITE),...
            SETTINGS.EMIN,...
            SETTINGS.EMAX,...
            SETTINGS.SYLVESTER_MATRIX_TYPE...
            );
        % 19 arguments
        
    end

    function WriteHeader()
        fprintf(fileID,'DATE,EX_NUM,m,n,o,t,ERROR_UXY,ERROR_VXY,ERROR_WXY,ERROR_DXY,MEAN_METHOD,BOOL_ALPHA_THETA,LOW_RANK_APPROX_METHOD,LRA_ITE,APF_METHOD,APF_ITE,error_min,error_max,Sylvester_Matrix_Type \n');
    end

end

function [dist_fxy] = GetError(fxy,fxy_exact)
% GetError : Get distance between f(x,y) and exact form.
%
% % Inputs.
%
% fxy : (Matrix) Coefficients of f(x,y) as computed.
%
% fxy_exact : (Matrix) Coefficients of f(x,y) exact.
%
% % Outputs.
%
% dist_fxy : Distance between f(x,y) and exact f(x,y)

dist_fxy = GetDistance(fxy_exact,fxy);

end

