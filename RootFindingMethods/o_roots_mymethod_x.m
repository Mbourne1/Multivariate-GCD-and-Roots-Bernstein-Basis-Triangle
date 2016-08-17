function [arr_wxy,vDegt_wx] = o_roots_mymethod_x(fxy_matrix,M)
% Given a bivariate polynomial compute the roots by differentiating with
% respect to x.

global SETTINGS

% Set the iteration number
ite = 1;

% Set the first entry of q to be the input polynomial f(x,y)
arr_fxy{1} = fxy_matrix;

% Get the total degree of f(x,y)
vDegt_arr_fxy(ite) = M(ite);

% Whilst the most recently calculated GCD has a degree greater than
% zero. ie is not a constant, perform a gcd calculation on it and its
% derivative.

while vDegt_arr_fxy(ite) > 0
    
    if (vDegt_arr_fxy(ite) == 1)
        % The derivative is a constant
        
        % The GCD is a constant
        arr_fxy{ite+1,1} = Differentiate_wrt_x(arr_fxy{ite,1});
        
        % Deconvolve 
        arr_uxy{ite+1,1} = Deconvolve_Bivariate(arr_fxy{ite,1},arr_fxy{ite+1,1});
        
        % Get total degree of d(x,y) and degree with respect to x and y
        vDegt_arr_fxy(ite+1,1) = vDegt_arr_fxy(ite,1)-1;
        
              
        break;
    end
    
    LineBreakLarge()
    fprintf([mfilename ' : ' sprintf('GCD Calculation Loop iteration = %i \n', ite)]);
    fprintf([mfilename ' : ' sprintf('Compute GCD of f_{%i} and derivative f_{%i}\n\n',ite,ite)]);
    
    
    
    % Get the total degree of f(x,y)
    m =  vDegt_arr_fxy(ite,1);
    
    % Get the derivative of f(x,y) with respect to x.
    gxy = Differentiate_wrt_x(arr_fxy{ite});
    
    % Get the total degree of g(x,y)
    n =  m - 1;
    
    % Get the upper and lower limit of the degree of the GCD(f, f')
    if ite > 1
        
        lower_lim = vDegt_arr_fxy(ite,1)- vNumDistinctRoots_fxy(ite-1,1);
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
    vDegt_arr_fxy(ite+1,1) = t;
    
    % Get number of distinct roots of f_{i}(x,y)
    vNumDistinctRoots_fxy(ite,1) = vDegt_arr_fxy(ite,1) - vDegt_arr_fxy(ite+1,1);
    
    fprintf([mfilename ' : ' sprintf('The computed deg(GCD(f_{%i},f_{%i}) is : %i \n',ite,ite,vDegt_arr_fxy(ite+1))]);
    fprintf([mfilename ' : ' sprintf('Number of distinct roots in f_{%i} : %i \n',ite,vNumDistinctRoots_fxy(ite))]);
    fprintf([mfilename ' : ' sprintf('Degree of f_{%i} : %i \n',ite + 1, vDegt_arr_fxy(ite+1))])
    
    LineBreakLarge()
    
    % Increment the iteration number
    ite = ite + 1;
    
end

% % Obtain the series h_{i}
% Each h_{x}(i) is obtained by the deconvolution of q_{x}(i) and q_{x}(i+1)

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
        
        arr_hxy = Get_arr_hxy(arr_fxy,vDegt_arr_fxy);
        
    case 'From ux'
        
        arr_hxy = arr_uxy;
        
    otherwise
        error('err')
        
end

% Get total degree of polynomials h(x,y)
vDegt_hxy = vDegt_arr_fxy(1:end-1) - vDegt_arr_fxy(2:end);


% %
% %
% Obtain the series of polynomials w_{x}{i}

% Each w_{x}(i) is obtained by the deconvolution of h_{x}(i) and h_{x}(i+1)

% Get number of polynomials in the array of h_{x}
nEntries_hxy = size(arr_hxy,1);

if nEntries_hxy > 1
    
    % Get array of polynomials w_{i}(x,y)
    arr_wxy = Get_arr_wxy(arr_hxy,vDegt_hxy);
  
    % Get total degree of polynomials w_{i}(x,y)
    vDegt_wx = vDegt_hxy(1:end-1) - vDegt_hxy(2:end);
    
    % Set the final w_{x}(i+1) to be equal to h_{x}(i+1)
    arr_wxy{end+1} = arr_hxy{end};
    
    % Set the final degree structure
    vDegt_wx(end) = vDegt_hxy(end);
    
else
    
    % Number of h_{i} is equal to one, so no deconvolutions to be performed
    arr_wxy{1} = arr_hxy{1};
    
    % Get the degree structure of h_{x,i}
    vDegt_wx(1) = vDegt_hxy(1);
end

for i = 1:1:length(arr_wxy)
    
    fprintf([mfilename ' : ' sprintf('Roots of degree %i',i) ' \n']);
    factor = arr_wxy{i};
    
    if (length(factor) > 1)
        
        display(factor./factor(2));
        
    end
    
end

end
