function [arr_wxy] = Get_arr_wxy(arr_hxy,vDegt_hxy)
%
%
% % Inputs
%
% arr_hxy : Array of polynomials h_{i}(x,y)
%
% vDegt_hxy : Vector of total degrees of h_{i}(x,y)

global SETTINGS

switch SETTINGS.WXY_DECONVOLUTION_METHOD
    
    case 'Separate' % Separate deconvolution
        
        arr_wxy = Deconvolve_Bivariate_Separate(arr_hxy);

        
    case 'Batch Without STLN'
        
        arr_wxy = Deconvolve_Bivariate_Batch(arr_hxy,vDegt_hxy);
        
    case 'Batch With STLN'
        
        arr_wxy = Deconvolve_Bivariate_Batch_With_STLN(arr_hxy,vDegt_hxy);
        
    otherwise
        
        display(SETTINGS.WXY_DECONVOLUTION_METHOD)
        
        Options_str = ...
            sprintf([
            '*Separate \n'...
            '*Batch Without STLN \n'...
            '*Batch With STLN \n'...
            '*Batch Constrained Without STLN \n'...
            '*Batch Constrained With STLN \n'
            ]);
        error([mfilename ' : ' sprintf('Deconvolution method must be one of the following : \n') Options_str])
        
end
end