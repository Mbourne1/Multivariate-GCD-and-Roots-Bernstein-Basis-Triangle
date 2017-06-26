function [arr_wxy] = Get_arr_wxy(arr_hxy,vDegt_hxy)
%
%
% % Inputs
%
% arr_hxy : Array of polynomials h_{i}(x,y)
%
% vDegt_hxy : Vector of total degrees of h_{i}(x,y)

global SETTINGS

switch SETTINGS.DECONVOLUTION_METHOD_WXY
    
    case 'Separate' % Separate deconvolution
        
        arr_wxy = Deconvolve_Bivariate_Separate(arr_hxy);

        
    case 'Batch Without STLN'
        
        arr_wxy = Deconvolve_Bivariate_Batch(arr_hxy,vDegt_hxy);
        
    case 'Batch With STLN'
        
        arr_wxy = Deconvolve_Bivariate_Batch_With_STLN(arr_hxy,vDegt_hxy);
        
    otherwise
        
        display(SETTINGS.DECONVOLUTION_METHOD_WXY)
        
        Options_str = ...
            sprintf([
            '*Separate \n'...
            '*Batch \n'...
            '*Batch With STLN \n'...
            '*Batch Constrained \n'...
            '*Batch Constrained With STLN \n'
            ]);
        error([mfilename ' : ' sprintf('Deconvolution method must be one of the following : \n') Options_str])
        
end
end