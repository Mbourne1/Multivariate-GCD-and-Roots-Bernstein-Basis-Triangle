function [t, GM_fx, GM_gx, alpha, th1, th2] = GetGCDDegree_2Polys(fxy, gxy, m, n, limits_t)
% GetGCDDegree(fxy,gxy,m,n)
%
% Inputs.
%
% [fxy, gxy] : Coefficients of polynomial f(x,y) and g(x,y)
%
% [m, n] : Total degree of f(x,y) and g(x,y)
%
% % Outputs
%
% t :
%
% [lambda, mu] :
%
% alpha :
%
% [th1, th2] :

global SETTINGS

my_limits = [0 min(m,n)];

lowerLimit = my_limits(1);
upperLimit = my_limits(2);

nSubresultants = upperLimit - lowerLimit + 1;


arr_SingularValues = cell(nSubresultants,1);
arr_R1 = cell(nSubresultants,1);

vGM_fx = zeros(nSubresultants,1);
vGM_gx = zeros(nSubresultants,1);

vAlpha = zeros(nSubresultants,1);
vTh1 = zeros(nSubresultants,1);
vTh2 = zeros(nSubresultants,1);

for i = 1:1:nSubresultants
    
    k = lowerLimit + (i-1);
    
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
        
        plotRowNorms(arr_R1RowNorms, my_limits, limits_t)
        plotMaxMinRowNorms(vMaxRowNormR1,vMinRowNormR1, my_limits, limits_t)
        
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
        
        
        plotDiagonalsR1(arr_R1, my_limits, limits_t)
        plotMaxMinDiagonalR1(vRatio_MaxMin_DiagonalEntry, my_limits, limits_t);
        
        
        metric = vRatio_MaxMin_DiagonalEntry;
        
    case 'Singular Values'
        
        vMinimumSingularValues = zeros(length(arr_SingularValues),1);
        for i = 1:1:length(arr_SingularValues)
           
            vec = arr_SingularValues{i};
            vMinimumSingularValues(i) = min(vec);
            
            
        end
        
        plotSingularValues(arr_SingularValues, my_limits, limits_t);
        plotMinimumSingularValues(vMinimumSingularValues, my_limits, limits_t);
        
        metric = vMinimumSingularValues;
        
    case 'Residuals'
        error('err')
    otherwise
        error('err')
end

if min(m,n) == 1
    
    t = GetGCDDegree_OneSubresultant(metric);
    
else
    
    t = GetGCDDegree_MultipleSubresultants(metric,[lowerLimit, upperLimit]);
    
    GM_fx = vGM_fx(t-lowerLimit + 1);
    GM_gx = vGM_gx(t-lowerLimit + 1);
    alpha = vAlpha(t-lowerLimit + 1);
    th1 = vTh1(t-lowerLimit + 1);
    th2 = vTh2(t-lowerLimit + 1);
    
end

end