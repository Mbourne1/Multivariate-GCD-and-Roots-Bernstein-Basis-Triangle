function [] = o_gcd(ex_num, emin, emax, mean_method, bool_alpha_theta, ...
    low_rank_approx_method, apf_method, sylvester_matrix_type)
% O_GCD : Compute the GCD of two polynomials f(x,y) and g(x,y) where f(x,y)
% and g(x,y) are given in the Bernstein form, and are taken from the
% example file.
%
% % Inputs.
%
% ex_num : Example Number
%
% el : Noise level - lower level
%
% em : Noise level - maximum level
%
% mean_method : Options for mean calculation method
%               *Geometric Mean Matlab Method
%               *Geometric Mean My Method
%               *None
%
% bool_alpha_theta : 'y'/'n'
%
% low_rank_approx_method : Options are
%               * Standard STLN
%               * None
%
% apf_method :
%       'None'
%       'Standard Linear APF'
%       'Standard Nonlinear APF'
%
% sylvester_matrix_type :
%       * T
%       * DT
%       * DTQ
%       * TQ
%
% % Example
%
% >> o_gcd('1',1e-10,1e-12,'Geometric Mean Matlab Method','y','Standard STLN','None','DTQ')


% % Ensure that minimum noise level is less than maximum noise level
if emin > emax
    temp = emin;
    emin = emax;
    emax = temp;
end

problem_type = 'GCD';

% % Set global variables
SetGlobalVariables(problem_type, ex_num, emin, emax, mean_method, bool_alpha_theta,...
    low_rank_approx_method, apf_method, sylvester_matrix_type);


% Add subfolders
restoredefaultpath

addpath (...
    'Bernstein Functions',...
    'Build Matrices',...
    'Build Sylvester Matrix',...
    'Formatting',...
    'Get Cofactors',...
    'Get GCD Coefficients',...
    'Get GCD Degree',...
    'Preprocessing'...
    );

addpath(genpath('APF'));
addpath(genpath('Examples'));
addpath(genpath('Low Rank Approximation'));

% Print parameters to console
fprintf('INPUTS. \n')
fprintf('EXAMPLE NUMBER %s \n',ex_num)
fprintf('EMIN : %s \n',emin)
fprintf('EMAX : %s \n',emax)
fprintf('MEAN METHOD : %s \n', mean_method)
fprintf('PREPROCESSING : %s \n',bool_alpha_theta)
fprintf('LOW RANK METHOD : %s \n',low_rank_approx_method)
fprintf('APF METHOD : %s \n', apf_method)


% %
% Get two polynomials f(x,y) and g(x,y) from example file.
[fxy,gxy,uxy_exact,vxy_exact,dxy_exact,m,n,t_exact] = Examples_GCD(ex_num);

% %
% Add noise to the coefficients

fxy = AddVariableNoiseToPoly(fxy,emin,emax);
gxy = AddVariableNoiseToPoly(gxy,emin,emax);

% %
% Get GCD by my method
[fxy,gxy,dxy,uxy, vxy, t] = o_gcd_mymethod(fxy,gxy,m,n);


% Get Error in u(x,y), v(x,y) and d(x,y)
dist_uxy = GetError(uxy,uxy_exact);
dist_vxy = GetError(vxy,vxy_exact);
dist_dxy = GetError(dxy,dxy_exact);

% Print the error in u(x,y), v(x,y) and d(x,y)
fprintf([mfilename ' : ' sprintf('Distance u(x,y) : %e \n', dist_uxy)]);
fprintf([mfilename ' : ' sprintf('Distance v(x,y) : %e \n', dist_vxy)]);
fprintf([mfilename ' : ' sprintf('Distance d(x,y) : %e \n', dist_dxy)]);

PrintToResultsFile(m,n,t,dist_dxy)

end

function [dist] = GetDistance(fxy,gxy)

fxy = fxy./fxy(1,1);
gxy = gxy./gxy(1,1);
try
    dist = norm(fxy - gxy) ./ norm(fxy);
catch
    dist = 1000;
end
end

function []= PrintToResultsFile(m,n,t,error_dx)

global SETTINGS

fullFileName = sprintf('Results/Results_o_gcd%s.txt',datetime('today'));

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
        
        % 15 FIELDS
        fprintf(fileID,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s \n',...
            datetime('now'),...
            SETTINGS.EX_NUM,...
            int2str(m),...
            int2str(n),...
            int2str(t),...
            error_dx,...
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
        
    end

    function WriteHeader()
        fprintf(fileID,'DATE,EX_NUM,m,n,t,ERROR_DXY,MEAN_METHOD,BOOL_ALPHA_THETA,LOW_RANK_APPROX_METHOD,LRA_ITE,APF_METHOD,APF_ITE,error_min,error_max,Sylvester_Matrix_Type \n');
    end

end

function [dist_fxy] = GetError(fxy,fxy_exact)
% GetError : Get distance between f(x,y) and exact form.
%
% % Inputs.
%
% fxy : Coefficients of f(x,y) as computed.
%
% fxy_exact : Coefficients of f(x,y) exact.
%
% % Outputs.
%
% dist_fxy : Distance between f(x,y) and exact f(x,y)

dist_fxy = GetDistance(fxy_exact,fxy);

end