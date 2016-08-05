function [] = o_gcd(ex_num,el,mean_method,bool_alpha_theta)
% Compute the GCD of two polynomials f(x,y) and g(x,y)
%
% % Inputs.
%
% ex_num : Example Number
%
% el : Noise level low
%
% mean_method : 
%               'Geometric Mean Matlab Method' 
%               'Geometric Mean My Method'
%               'None'
%
% bool_alpha_theta : 'y'/'n' 
% 
% % Example
%
% >> o_gcd('1',1e-10,'Geometric Mean Matlab Method','y')


addpath (...
        'Examples',...
        'GetGCDDegree',...
        'Preprocessing'...
        );
    
global SETTINGS
SETTINGS.EX_NUM = ex_num;
SETTINGS.EMIN = el;

SetGlobalVariables(mean_method,bool_alpha_theta);

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


display(uxy_exact)
display(uxy)
display(uxy_exact./uxy);

display(vxy_exact)
display(vxy)
display(vxy_exact./vxy);


display(dxy_exact)
display(dxy)
display(dxy_exact./dxy);

dist_dxy = GetDistance(dxy_exact,dxy);
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