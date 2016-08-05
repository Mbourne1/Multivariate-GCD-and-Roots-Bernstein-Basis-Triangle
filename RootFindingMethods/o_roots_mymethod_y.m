function [wy,vDegt_wy] = o_roots_mymethod_y(fxy_matrix,M)
% Given a bivariate polynomial compute the roots by differentiating with
% respect to y.

global SETTINGS

% Set the iteration number
ite = 1;

% Set the first entry of q to be the input polynomial f(x,y)
fy{1} = fxy_matrix;

% Get the total degree of f(x,y)
vDegt_fy(ite) = M(ite);

% Whilst the most recently calculated GCD has a degree greater than
% zero. ie is not a constant, perform a gcd calculation on it and its
% derivative.
while vDegt_fy(ite) > 0
    
    if (vDegt_fy(ite) == 1)
        % The derivative is a constant
        
        % The GCD is a constant
        fy{ite+1} = Differentiate_wrt_y(fy{ite});
        
        % Deconvolve
        uy{ite+1} = Deconvolve_Bivariate_Single(fy{ite},fy{ite+1});
        
        % Get total degree of d(x,y) and degree with respect to x and y
        vDegt_fy(ite+1) = vDegt_fy(ite)-1;
        
        break;
    end
    
    LineBreakLarge()
    fprintf([mfilename ' : ' sprintf('GCD calculation wrt y iteration : %i \n\n',ite)]);
    fprintf([mfilename ' : ' sprintf('Compute GCD of f_{%i} and derivative f_{%i}\n\n',ite,ite)]);
    
    % Get the derivative of f(x,y) with respect to y.
    gxy = Differentiate_wrt_y(fy{ite});
    
    % Get the total degree of f(x,y)
    m = vDegt_fy(ite);
    
    % Get the total degree of g(x,y)
    n = m - 1;
    
    if ite > 1
        lower_lim = vDegt_fy(ite)-dy(ite-1);
        upper_lim = m-1;
    else
        lower_lim = 1;
        upper_lim = m-1;
    end
    
    
    fprintf([mfilename ' : ' sprintf('Minimum degree of f_{%i}: %i \n', ite+1, lower_lim)]);
    fprintf([mfilename ' : ' sprintf('Maximum degree of f_{%i}: %i \n', ite+1, upper_lim)]);
    
    % Get the GCD of f(x,y) and g(x,y)
    [fy{ite},~,fy{ite+1},uy{ite},~,t] = o_gcd_mymethod(fy{ite},gxy,m,n,[lower_lim,upper_lim]);
    
    % Get total degree of d(x,y) and degree with respect to x and y
    vDegt_fy(ite+1) = t;
    
    % Get number of distinct roots of f(ite)
    dy(ite) = vDegt_fy(ite) - vDegt_fy(ite+1);
    
    fprintf([mfilename ' : ' sprintf('The computed deg(GCD(f_{%i},f_{%i}) is : %i \n',ite,ite,vDegt_fy(ite+1))])
    fprintf([mfilename ' : ' sprintf('Number of distinct roots in f_{%i} : %i \n',ite,dy(ite))])
    fprintf([mfilename ' : ' sprintf('Degree of f_{%i} : %i \n',ite + 1, vDegt_fy(ite+1))])
    
    LineBreakLarge()
    
    % Increment the iteration number
    ite = ite + 1;
end

% %
% %
% % Obtain the series h_{i}


% Get number of elements in the series of polynomials q_{i}
[~,num_entries_fy] = size(fy);

% %
% %     Get h_{i}(y)
% %

% METHOD refers to method used to compute h{i}(x,y), either by
% deconvolution or from GCD triples computed above.
% From Deconvolutions
% From uy
SETTINGS.HXY_METHOD = 'From Deconvolutions';
switch SETTINGS.HXY_METHOD
    
    case 'From Deconvolutions'
        
        switch SETTINGS.DECONVOLUTION_METHOD
            
            case 'Separate'
                % For each pair of q_{x}(i) and q_{x}(i+1)
                for i = 1 : 1 : num_entries_fy - 1
                    
                    % Perform deconvolutions
                    hy{i} = Deconvolve_Bivariate_Single(fy{i},fy{i+1});
                    
                end
            case 'Batch'
                
                % Get the set of polynomials hy{i} from the deconvolution of the
                % set of polynomials fy{i}/fy{i+1}
                hy = Deconvolve_Bivariate_Batch(fy,vDegt_fy);
                 
            otherwise
                error('err')
                
        end
        
    case 'From ux'
        
        hx = ux;
        
    otherwise
        error('err');
        
end


vDegt_hy = vDegt_fy(1:end-1) - vDegt_fy(2:end);


% %
% %
% Obtain the series of polynomials w_{y}{i} for the

% Each w_{y}(i) is obtained by the deconvolution of h_{y}(i) and h_{y}(i+1)

% Get number of polynomials in the array of h_{y}
[~,num_entries_hy] = size(hy);

if num_entries_hy > 1
    
    switch SETTINGS.DECONVOLUTION_METHOD
        
        case 'Separate' % Separate deconvolution
            
            for i = 1 : 1 : num_entries_hy - 1 % For each pair of h_{y}(i) and h_{y}(i+1)
                
                % Deconvolve
                wy{i} = Deconvolve_Bivariate_Single(hy{i},hy{i+1});
                
            end

        case 'Batch' % Batch deconvolution
            wy = Deconvolve_Bivariate_Batch(hy,vDegt_hy);
        otherwise
            error('err')
            
    end
    

    vDegt_wy = vDegt_hy(1:end-1) - vDegt_hy(2:end);
    
    
    % Set the final w_{x}(i+1) to be equal to h_{x}(i+1)
    wy{end+1} = hy{end};
    
    % Set the final degree structure

    vDegt_wy(end) = vDegt_hy(end);
    
else
    
    % Number of h_{i} is equal to one, so no deconvolutions to be performed
    wy{1} = hy{1};
    
    % Get the degree structure of h_{x,i}

    vDegt_wy(1) = vDegt_hy(1);
end

for i = 1:1:length(wy)
    fprintf([mfilename ' : ' sprintf('Roots of degree %i',i) ' \n']);
    factor = wy{i};
    if (length(factor) > 1)
        display(factor./factor(2));
    end
    
end

end