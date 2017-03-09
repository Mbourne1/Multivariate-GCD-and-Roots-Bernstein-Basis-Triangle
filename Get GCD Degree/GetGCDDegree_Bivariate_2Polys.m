function [t, GM_fx, GM_gx, alpha, th1, th2] = GetGCDDegree_Bivariate_2Polys(fxy, gxy, m, n, limits_t)
% GetGCDDegree(fxy,gxy,m,n)
%
% Inputs.
%
% fxy : (Matrix) Coefficients of the polynomial f(x,y)
% gxy : (Matrix) Coefficients of the polynomial g(x,y)
% m : (Int) Total degree of f(x,y)
% n : (Int) Total degree of g(x,y)
%
% limits_t = []
%
% % Outputs
%
% t : (Int)
% GM_fx : (Float)
% GM_gx : (Float)
% alpha : (Float)
% th1 : (Float)
% th2 : (Float)

global SETTINGS

myLimits = [0 min(m,n)];

myLowerLimit = myLimits(1);
myUpperLimit = myLimits(2);

nSubresultants = myUpperLimit - myLowerLimit + 1;


arr_SingularValues = cell(nSubresultants,1);
arr_R1 = cell(nSubresultants,1);

vGM_fx = zeros(nSubresultants,1);
vGM_gx = zeros(nSubresultants,1);

vAlpha = zeros(nSubresultants,1);
vTh1 = zeros(nSubresultants,1);
vTh2 = zeros(nSubresultants,1);

for i = 1:1:nSubresultants
    
    k = myLowerLimit + (i-1);
    
    % %
    % Geometric Mean
    
    % Get Geometric mean of f(x,y) in C_{n-k}(f)
    vGM_fx(i) = GetGeometricMean(fxy, m, n-k);
    
    % Get Geometric mean of g(x,y) in C_{m-k}(g)
    vGM_gx(i) = GetGeometricMean(gxy, n, m-k);
    
    % Divide entries of f(x,y) and g(x,y) by geometric mean
    fxy_n = fxy ./ vGM_fx(i);
    gxy_n = gxy ./ vGM_gx(i);
    
    % %
    % Preprocess
    [vAlpha(i), vTh1(i), vTh2(i)] = Preprocess(fxy_n, gxy_n, m, n, k);
    
    % Get f(x,y) and g(x,y) with thetas to get f(w_{1},w_{2}) and g(w_{1},w_{2})
    fww = GetWithThetas(fxy_n, m, vTh1(i), vTh2(i));
    gww = GetWithThetas(gxy_n, n, vTh1(i), vTh2(i));
    
    % %
    % Build the Sylvester Matrix
    Sk = BuildSylvesterMatrix_2Polys(fww, vAlpha(i).*gww, m, n, k);
    
    % Get vector of singular values
    arr_SingularValues{i} = svd(Sk);
    
    [~,R] = qr(Sk);
    [~,c] = size(R);
    
    arr_R1{i} = R(1:c, 1:c);
    
end

% Metric used to compute the degree of the GCD
% R1 Row Norms
% R1 Row Diagonals
% Singular Values
% Residuals
switch SETTINGS.RANK_REVEALING_METRIC
    case 'R1 Row Norms'
        
        arr_R1_RowNorms = cell(nSubresultants,1);
        vMaxRowNormR1 = zeros(nSubresultants,1);
        vMinRowNormR1 = zeros(nSubresultants,1);
        
        for i = 1:1:nSubresultants
            
            arr_R1_RowNorms{i} = sqrt(sum(arr_R1{i}.^2,2))./norm(arr_R1{i});
            
            % Get maximum and minimum row norms of rows of R1.
            vMaxRowNormR1(i) = max(arr_R1_RowNorms{i});
            vMinRowNormR1(i) = min(arr_R1_RowNorms{i});
            
        end
        
        if(SETTINGS.PLOT_GRAPHS)
            plotRowNorms(arr_R1RowNorms, myLimits, limits_t)
            plotMaxMinRowNorms(vMaxRowNormR1, vMinRowNormR1, myLimits, limits_t)
        end
        metric = vMinRowNormR1 ./ vMaxRowNormR1;
        
    case 'R1 Row Diagonals'
        
        vMaxDiagonalEntry = zeros(nSubresultants,1);
        vMinDiagonalEntry = zeros(nSubresultants,1);
        
        for i = 1:1:length(arr_R1)
            
            % Get maximum diagonal
            vMaxDiagonalEntry(i) = max(abs(diag(arr_R1{i})));
            
            % Get minimum diagonal
            vMinDiagonalEntry(i) = min(abs(diag(arr_R1{i})));
            
        end
        
        vRatio_MaxMin_DiagonalEntry = vMinDiagonalEntry ./ vMaxDiagonalEntry;
        
        if(SETTINGS.PLOT_GRAPHS)
            %plotDiagonalsR1(arr_R1, my_limits, limits_t)
            plotMaxMinDiagonalR1(vRatio_MaxMin_DiagonalEntry, myLimits, limits_t);
        end
        
        metric = vRatio_MaxMin_DiagonalEntry;
        
    case 'Singular Values'
        
        vMinimumSingularValues = zeros(length(arr_SingularValues),1);
        for i = 1:1:length(arr_SingularValues)
            
            vec = arr_SingularValues{i};
            vMinimumSingularValues(i) = min(vec);
            
            
        end
        
        if(SETTINGS.PLOT_GRAPHS)
            plotSingularValues(arr_SingularValues, myLimits, limits_t);
            plotMinimumSingularValues(vMinimumSingularValues, myLimits, limits_t);
        end
        
        metric = vMinimumSingularValues;
        

    case 'Residuals'
        error('err')

    otherwise
        error('err')
end

if (myLowerLimit == myUpperLimit)
    
    t = GetGCDDegree_OneSubresultant(metric);
    
else
    
    t = GetGCDDegree_MultipleSubresultants(metric,[myLowerLimit, myUpperLimit]);
    
    GM_fx = vGM_fx(t-myLowerLimit + 1);
    GM_gx = vGM_gx(t-myLowerLimit + 1);
    alpha = vAlpha(t-myLowerLimit + 1);
    th1 = vTh1(t-myLowerLimit + 1);
    th2 = vTh2(t-myLowerLimit + 1);
    
end

end