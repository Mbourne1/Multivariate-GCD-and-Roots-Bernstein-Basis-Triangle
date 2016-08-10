function [wx,vDegt_wx] = o_roots_mymethod_x(fxy_matrix,M)
% Given a bivariate polynomial compute the roots by differentiating with
% respect to x.

global SETTINGS

% Set the iteration number
ite = 1;

% Set the first entry of q to be the input polynomial f(x,y)
arr_fxy{1} = fxy_matrix;

% Get the total degree of f(x,y)
vDegt_fxy(ite) = M(ite);

% Whilst the most recently calculated GCD has a degree greater than
% zero. ie is not a constant, perform a gcd calculation on it and its
% derivative.

while vDegt_fxy(ite) > 0
    
    if (vDegt_fxy(ite) == 1)
        % The derivative is a constant
        
        % The GCD is a constant
        arr_fxy{ite+1,1} = Differentiate_wrt_x(arr_fxy{ite,1});
        
        % Deconvolve 
        arr_uxy{ite+1,1} = Deconvolve_Bivariate(arr_fxy{ite,1},arr_fxy{ite+1,1});
        
        % Get total degree of d(x,y) and degree with respect to x and y
        vDegt_fxy(ite+1,1) = vDegt_fxy(ite,1)-1;
        
              
        break;
    end
    
    LineBreakLarge()
    fprintf([mfilename ' : ' sprintf('GCD Calculation Loop iteration = %i \n', ite)]);
    fprintf([mfilename ' : ' sprintf('Compute GCD of f_{%i} and derivative f_{%i}\n\n',ite,ite)]);
    
    
    
    % Get the total degree of f(x,y)
    m =  vDegt_fxy(ite,1);
    
    % Get the derivative of f(x,y) with respect to x.
    gxy = Differentiate_wrt_x(arr_fxy{ite});
    
    % Get the total degree of g(x,y)
    n =  m - 1;
    
    % Get the upper and lower limit of the degree of the GCD(f, f')
    if ite > 1
        
        lower_lim = vDegt_fxy(ite,1)- vNumDistinctRoots_fxy(ite-1,1);
        upper_lim = m-1;
    else
        lower_lim = 1;
        upper_lim = m-1;
    end
    
    fprintf([mfilename ' : ' sprintf('Minimum degree of f_{%i}: %i \n', ite+1, lower_lim)]);
    fprintf([mfilename ' : ' sprintf('Maximum degree of f_{%i}: %i \n', ite+1, upper_lim)]);
    
    % GCD is only a scalar with respect to x so set equal to g(x,y).
    [arr_fxy{ite,1},~,arr_fxy{ite+1,1},arr_uxy{ite,1},vxy{ite,1},t] = o_gcd_mymethod(arr_fxy{ite,1},gxy,m,n,[lower_lim, upper_lim]);
    
    % Get total degree of d(x,y) and degree with respect to x and y
    vDegt_fxy(ite+1,1) = t;
    
    % Get number of distinct roots of f_{i}(x,y)
    vNumDistinctRoots_fxy(ite,1) = vDegt_fxy(ite,1) - vDegt_fxy(ite+1,1);
    
    fprintf([mfilename ' : ' sprintf('The computed deg(GCD(f_{%i},f_{%i}) is : %i \n',ite,ite,vDegt_fxy(ite+1))]);
    fprintf([mfilename ' : ' sprintf('Number of distinct roots in f_{%i} : %i \n',ite,vNumDistinctRoots_fxy(ite))]);
    fprintf([mfilename ' : ' sprintf('Degree of f_{%i} : %i \n',ite + 1, vDegt_fxy(ite+1))])
    
    LineBreakLarge()
    
    % Increment the iteration number
    ite = ite + 1;
    
end

% % Obtain the series h_{i}
% Each h_{x}(i) is obtained by the deconvolution of q_{x}(i) and q_{x}(i+1)

% Get number of polynomials in the series of polynomials q_{i}
nEntries_arr_fxy = size(arr_fxy,1);

% %
% %     Get h_{i}(x,y)
% %

% METHOD refers to method used to compute h{i}(x,y), either by
% deconvolution or from GCD triples computed above.
% From Deconvolutions
% From ux

SETTINGS.HXY_METHOD = 'From Deconvolutions';
switch SETTINGS.HXY_METHOD
    
    case 'From Deconvolutions'
        
        switch SETTINGS.DECONVOLUTION_METHOD
            
            case 'Separate' % Separate deconvolution
                
                arr_hxy = cell(nEntries_arr_fxy - 1,1);
                
                for i = 1 : 1 : nEntries_arr_fxy - 1 % For each pair of q_{x}(i) and q_{x}(i+1)
                    
                    % Deconvolve
                    arr_hxy{i,1} = Deconvolve_Bivariate_Single(arr_fxy{i}, arr_fxy{i+1});
                    
                end
            case 'Batch' % Batch deconvolution
                
                % Get the set of polynomials hx{i} from the deconvolution of the
                % set of polynomials fx{i}/fx{i+1}
                arr_hxy = Deconvolve_Bivariate_Batch(arr_fxy,vDegt_fxy);
                
            case 'Batch Constrained'
                
                % Get the set of polynomials hx{i} from the deconvolution of the
                % set of polynomials fx{i}/fx{i+1}
                arr_hxy = Deconvolve_Bivariate_Batch_Constrained(arr_fxy,vDegt_fxy);
                
                
            otherwise
                error([mfilename ' : ' sprintf(' Deconvolution Method is either (Separate) or (Batch) or (Batch Constrained)')])
                
        end
        
    case 'From ux'
        
        arr_hxy = arr_uxy;
        
    otherwise
        error('err')
        
end

vDegt_hx = vDegt_fxy(1:end-1) - vDegt_fxy(2:end);


% %
% %
% Obtain the series of polynomials w_{x}{i}

% Each w_{x}(i) is obtained by the deconvolution of h_{x}(i) and h_{x}(i+1)

% Get number of polynomials in the array of h_{x}
nEntries_hxy = size(arr_hxy,1);

if nEntries_hxy > 1
    
    switch SETTINGS.DECONVOLUTION_METHOD
        
        case 'Separate' % Separate deconvolution
            
            for i = 1 : 1 : nEntries_hxy - 1 % For each pair of q_{x}(i) and q_{x}(i+1)
                
                % Deconvolve
                wx{i} = Deconvolve_Bivariate_Single(arr_hxy{i},arr_hxy{i+1});
                
            end

        case 'Batch' % Batch deconvolution
            wx = Deconvolve_Bivariate_Batch(arr_hxy,vDegt_hx);
            
        case 'Batch Constrained'
            wx = Deconvolve_Bivariate_Batch(arr_hxy,vDegt_hx);
        otherwise
            error('err')
            
    end
    
    
    vDegt_wx = vDegt_hx(1:end-1) - vDegt_hx(2:end);
    
    
    % Set the final w_{x}(i+1) to be equal to h_{x}(i+1)
    wx{end+1} = arr_hxy{end};
    
    % Set the final degree structure
    vDegt_wx(end) = vDegt_hx(end);
    
else
    
    % Number of h_{i} is equal to one, so no deconvolutions to be performed
    wx{1} = arr_hxy{1};
    % Get the degree structure of h_{x,i}
    vDegt_wx(1) = vDegt_hx(1);
end

for i = 1:1:length(wx)
    fprintf([mfilename ' : ' sprintf('Roots of degree %i',i) ' \n']);
    factor = wx{i};
    
    if (length(factor) > 1)
        
        display(factor./factor(2));
        
    end
    
end

end
