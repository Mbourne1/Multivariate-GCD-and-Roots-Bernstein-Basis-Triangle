function [wx,vDegt_wx] = o_roots_mymethod_x(fxy_matrix,M)
% Given a bivariate polynomial compute the roots by differentiating with
% respect to x.

global SETTINGS

% Set the iteration number
ite = 1;

% Set the first entry of q to be the input polynomial f(x,y)
fx{1} = fxy_matrix;

% Get the total degree of f(x,y)
vDegt_fx(ite) = M(ite);

% Whilst the most recently calculated GCD has a degree greater than
% zero. ie is not a constant, perform a gcd calculation on it and its
% derivative.

while vDegt_fx(ite) > 0
    
    if (vDegt_fx(ite) == 1)
        % The derivative is a constant
        
        % The GCD is a constant
        fx{ite+1} = Differentiate_wrt_x(fx{ite});
        
        % Deconvolve
        ux{ite+1} = Deconvolve_Bivariate_Single(fx{ite},fx{ite+1});
        
        % Get total degree of d(x,y) and degree with respect to x and y
        vDegt_fx(ite+1) = vDegt_fx(ite)-1;
        
        break;
    end
    
    LineBreakLarge()
    fprintf([mfilename ' : ' sprintf('GCD Calculation Loop iteration = %i \n', ite)]);
    fprintf([mfilename ' : ' sprintf('Compute GCD of f_{%i} and derivative f_{%i}\n\n',ite,ite)]);
    
    
    
    % Get the total degree of f(x,y)
    m =  vDegt_fx(ite);
    
    % Get the derivative of f(x,y) with respect to x.
    gxy = Differentiate_wrt_x(fx{ite});
    
    % Get the total degree of g(x,y)
    n =  m - 1;
    
    % Get the upper and lower limit of the degree of the GCD(f, f')
    if ite > 1
        lower_lim = vDegt_fx(ite)-dx(ite-1);
        upper_lim = m-1;
    else
        lower_lim = 1;
        upper_lim = m-1;
    end
    
    fprintf([mfilename ' : ' sprintf('Minimum degree of f_{%i}: %i \n', ite+1, lower_lim)]);
    fprintf([mfilename ' : ' sprintf('Maximum degree of f_{%i}: %i \n', ite+1, upper_lim)]);
    
    % GCD is only a scalar with respect to x so set equal to g(x,y).
    [fx{ite},~,fx{ite+1},ux{ite},vx{ite},t] = o_gcd_mymethod(fx{ite},gxy,m,n,[lower_lim, upper_lim]);
    
    % Get total degree of d(x,y) and degree with respect to x and y
    vDegt_fx(ite+1) = t;
    
    % Get number of distinct roots of f(ite)
    dx(ite) = vDegt_fx(ite) - vDegt_fx(ite+1);
    
    fprintf([mfilename ' : ' sprintf('The computed deg(GCD(f_{%i},f_{%i}) is : %i \n',ite,ite,vDegt_fx(ite+1))]);
    fprintf([mfilename ' : ' sprintf('Number of distinct roots in f_{%i} : %i \n',ite,dx(ite))]);
    fprintf([mfilename ' : ' sprintf('Degree of f_{%i} : %i \n',ite + 1, vDegt_fx(ite+1))])
    
    LineBreakLarge()
    
    % Increment the iteration number
    ite = ite + 1;
    
end

% % Obtain the series h_{i}
% Each h_{x}(i) is obtained by the deconvolution of q_{x}(i) and q_{x}(i+1)

% Get number of polynomials in the series of polynomials q_{i}
[~,nEntries_fx] = size(fx);

% %
% %     Get h_{i}(x)
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
                
                for i = 1 : 1 : nEntries_fx - 1 % For each pair of q_{x}(i) and q_{x}(i+1)
                    
                    % Deconvolve
                    hx{i} = Deconvolve_Bivariate_Single(fx{i}, fx{i+1});
                    
                end
            case 'Batch' % Batch deconvolution
                
                % Get the set of polynomials hx{i} from the deconvolution of the
                % set of polynomials fx{i}/fx{i+1}
                hx = Deconvolve_Bivariate_Batch(fx,vDegt_fx);
                
            otherwise
                error([mfilename ' : ' sprintf(' Deconvolution Method is either Separate or Batch')])
                
        end
        
    case 'From ux'
        
        hx = ux;
        
    otherwise
        error('err')
        
end

vDegt_hx = vDegt_fx(1:end-1) - vDegt_fx(2:end);


% %
% %
% Obtain the series of polynomials w_{x}{i}

% Each w_{x}(i) is obtained by the deconvolution of h_{x}(i) and h_{x}(i+1)

% Get number of polynomials in the array of h_{x}
[~,nEntries_hx] = size(hx);

if nEntries_hx > 1
    
    switch SETTINGS.DECONVOLUTION_METHOD
        
        case 'Separate' % Separate deconvolution
            
            for i = 1 : 1 : nEntries_hx - 1 % For each pair of q_{x}(i) and q_{x}(i+1)
                
                % Deconvolve
                wx{i} = Deconvolve_Bivariate_Single(hx{i},hx{i+1});
                
            end

        case 'Batch' % Batch deconvolution
            wx = Deconvolve_Bivariate_Batch(hx,vDegt_hx);
        otherwise
            error('err')
            
    end
    
    
    vDegt_wx = vDegt_hx(1:end-1) - vDegt_hx(2:end);
    
    
    % Set the final w_{x}(i+1) to be equal to h_{x}(i+1)
    wx{end+1} = hx{end};
    
    % Set the final degree structure
    vDegt_wx(end) = vDegt_hx(end);
    
else
    
    % Number of h_{i} is equal to one, so no deconvolutions to be performed
    wx{1} = hx{1};
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
