function [] = o_gcd(ex_num, el, em, mean_method, bool_alpha_theta, ...
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


% Set variables
global SETTINGS
SETTINGS.PLOT_GRAPHS = 'y';
SETTINGS.EX_NUM = ex_num;

if el > em
    temp = el;
    el = em;
    em = temp;
end

SETTINGS.EMIN = el;
SETTINGS.EMAX = em;




% Add subfolders
restoredefaultpath

addpath (...
    'Bernstein Functions',...
    'Build Matrices',...
    'Formatting',...
    'Get Cofactors',...
    'Get GCD Coefficients',...
    'Get GCD Degree',...
    'Preprocessing'...
    );
addpath(genpath('APF'));
addpath(genpath('Examples'));
addpath(genpath('Low Rank Approximation'));



% % Set global variables
SetGlobalVariables(mean_method, bool_alpha_theta,...
    low_rank_approx_method, apf_method, sylvester_matrix_type);

% %
% Get two polynomials f(x,y) and g(x,y) from example file.
[fxy,gxy,uxy_exact,vxy_exact,dxy_exact,m,n,t_exact] = Examples_GCD(ex_num);

% %
% Add noise to the coefficients

fxy = AddVariableNoiseToPoly(fxy,el,em);
gxy = AddVariableNoiseToPoly(gxy,el,em);

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

fullFileName = 'Results/Results_GCD.txt';


if exist(fullFileName, 'file')
    fileID = fopen(fullFileName,'a');
    fprintf(fileID,'%s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s \n',...
        datetime('now'),...
        SETTINGS.EX_NUM,...
        int2str(m),...
        int2str(n),...
        int2str(t),...
        error_dx,...
        SETTINGS.MEAN_METHOD,...
        SETTINGS.BOOL_ALPHA_THETA,...
        SETTINGS.LOW_RANK_APPROX_METHOD,...
        SETTINGS.APF_METHOD,...
        SETTINGS.EMIN,...
        SETTINGS.EMAX,...
        SETTINGS.SYLVESTER_MATRIX_TYPE...
        );
    fclose(fileID);
else
    % File does not exist.
    warningMessage = sprintf('Warning: file does not exist:\n%s', fullFileName);
    uiwait(msgbox(warningMessage));
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