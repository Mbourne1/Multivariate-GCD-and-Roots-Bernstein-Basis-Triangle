function arr_hxy = Get_arr_hxy(arr_fxy,vDegt_fxy)
% GET_ARR_HXY
%
% Given the set of polynomials f(x,y) compute h(x,y).
%
% % Inputs
%
% arr_fxy : Array of polynomials f(x,y)
%
% vDegt_fxy : Vector of Total degree of polynomials f(x,y)

global SETTINGS

switch SETTINGS.HXY_DECONVOLUTION_METHOD
    
    case 'Separate' % Separate deconvolution
        
        arr_hxy = Deconvolve_Bivariate_Separate(arr_fxy);
        
    case 'Batch Without STLN' % Batch deconvolution
        
        % Get the set of polynomials hx{i} from the deconvolution of the
        % set of polynomials fx{i}/fx{i+1}
        arr_hxy = Deconvolve_Bivariate_Batch_Without_STLN(arr_fxy,vDegt_fxy);
        
    case 'Batch With STLN'
        
        arr_hxy = Deconvolve_Bivariate_Batch_With_STLN(arr_fxy,vDegt_fxy);
        
    case 'Batch Constrained Without STLN'
        
        % Get the set of polynomials hx{i} from the deconvolution of the
        % set of polynomials fx{i}/fx{i+1}
        arr_hxy = Deconvolve_Bivariate_Batch_Constrained_Without_STLN(arr_fxy,vDegt_fxy);
        
    case 'Batch Constrained With STLN'
        
        arr_hxy = Deconvolve_Bivariate_Batch_Constrained_With_STLN(arr_fxy,vDegt_fxy);
        
    otherwise
        
        display(SETTINGS.HXY_DECONVOLUTION_METHOD)
        
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