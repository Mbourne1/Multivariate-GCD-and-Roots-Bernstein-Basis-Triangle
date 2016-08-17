function [] = o_gcd(ex_num,el,mean_method,bool_alpha_theta,low_rank_approx_method)
% O_GCD : Compute the GCD of two polynomials f(x,y) and g(x,y)
%
% % Inputs.
%
% ex_num : Example Number
%
% el : Noise level low
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
%
% % Example
%
% >> o_gcd('1',1e-10,'Geometric Mean Matlab Method','y','Standard STLN')


addpath (...
        'Examples',...
        'Formatting',...
        'GetCofactors',...
        'GetGCDCoefficients',...
        'GetGCDDegree',...
        'Low Rank Approximation',...
        'Preprocessing'...
        );
    
global SETTINGS
SETTINGS.EX_NUM = ex_num;
SETTINGS.EMIN = el;

SetGlobalVariables(mean_method,bool_alpha_theta,low_rank_approx_method);

% %
% Get two polynomials f(x,y) and g(x,y) from example file.
[fxy,gxy,uxy_exact,vxy_exact,dxy_exact,m,n,t_exact] = Examples_GCD(ex_num);

% %
% Add noise to the coefficients

fxy = Noise(fxy,el);
gxy = Noise(gxy,el);

% %
% Get GCD by my method
[fxy,gxy,dxy,uxy, vxy, t] = o_gcd_mymethod(fxy,gxy,m,n);


dist_uxy = GetError(uxy,uxy_exact);
dist_vxy = GetError(vxy,vxy_exact);
dist_dxy = GetError(dxy,dxy_exact);

fprintf([mfilename ' : ' sprintf('Distance u(x,y) : %e \n', dist_uxy)]);
fprintf([mfilename ' : ' sprintf('Distance v(x,y) : %e \n', dist_vxy)]);
fprintf([mfilename ' : ' sprintf('Distance d(x,y) : %e \n', dist_dxy)]);

PrintToResultsFile(m,n,t,dist_dxy)

end

function [dist] = GetDistance(fxy,gxy)

fxy = fxy./fxy(1,1);
gxy = gxy./gxy(1,1);

dist = norm(fxy - gxy) ./ norm(fxy);

end

function []= PrintToResultsFile(m,n,t,error_dx)

global SETTINGS

fullFileName = 'Results_GCD.txt';


if exist(fullFileName, 'file')
    fileID = fopen(fullFileName,'a');
    fprintf(fileID,'%s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \n',...
        datetime('now'),...
        SETTINGS.EX_NUM,...
        int2str(m),...
        int2str(n),...
        int2str(t),...
        error_dx,...
        SETTINGS.MEAN_METHOD,...
        SETTINGS.BOOL_ALPHA_THETA,...
        SETTINGS.EMIN...
        );
    fclose(fileID);
else
    % File does not exist.
    warningMessage = sprintf('Warning: file does not exist:\n%s', fullFileName);
    uiwait(msgbox(warningMessage));
end
end

function [dist_fxy] = GetError(fxy,fxy_exact)
% GETERROR : Get distance between f(x,y) and exact form.
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