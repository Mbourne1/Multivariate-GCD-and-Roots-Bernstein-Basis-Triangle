function [t, GM_fx, GM_gx, GM_hx, alpha, th1, th2, rank_range] = ...
    GetGCDDegree_Bivariate_3Polys(fxy, gxy, hxy, m, n, o, limits_t, rank_range)
% GetGCDDegree_Bivariate_3Polys(fxy, gxy, hxy, m, n, o)
%
% % Inputs.
%
% fxy : (Vector) Coefficients of polynomial f(x,y)
%
% gxy : (Vector) Coefficients of polynomial g(x,y)
%
% hxy : (Vector) Coefficients of polynomial h(x,y)
%
% m : (Int) Degree of polynomial f(x,y)
%
% n : (Int) Degree of polynomial g(x,y)
%
% o : (Int) Degree of polynomial h(x,y)
%
% limits_t : [(Int) (Int)] 
%
% rank_range : [(Float) (Float) ]
%
% % Outputs
%
% t : (Int) Degree of GCD d(x,y)
%
% GM_fx : (Float) Geometric mean of coefficients of f(x)
%
% GM_gx : (Float) Geometric mean of coefficients of g(x)
%
% GM_hx : (Float) Geometric mean of coefficients of h(x)
%
% alpha : (Float) Optimal value of \alpha
%
% th1 : (Float) Optimal value of \theta_{1}
%
% th2 : (Float) Optimal value of \theta_{2}

global SETTINGS

% Get range of Sylvester subresultant matrices to be constructed
lowerLimit_k = 0;
upperLimit_k = min([m,n,o]);
limits_k = [lowerLimit_k upperLimit_k];

rank_range_low = rank_range(1);
rank_range_high = rank_range(2);

% Get number of Sylvester Subresultants
nSubresultants = upperLimit_k - lowerLimit_k + 1;

arr_SingularValues = cell(nSubresultants, 1);
arr_R1 = cell(nSubresultants, 1);
vGM_fx = zeros(nSubresultants, 1);
vGM_gx = zeros(nSubresultants, 1);
vGM_hx = zeros(nSubresultants, 1);
vAlpha = zeros(nSubresultants, 1);
vBeta = zeros(nSubresultants, 1);
vTh1 = zeros(nSubresultants, 1);
vTh2 = zeros(nSubresultants, 1);


for i = 1:1:nSubresultants
    
    k = lowerLimit_k + (i-1);
    
    % %
    % Geometric Mean
    
    % Get Geometric mean of f(x,y) in C_{n-k}(f), g(x,y) and h(x,y)
    vGM_fx(i) = GetMean(fxy, m, n-k);
    vGM_gx(i) = GetMean(gxy, n, m-k);
    vGM_hx(i) = GetMean(hxy, o, n-k);
    
    % Divide entries of f(x,y) and g(x,y) by geometric mean
    fxy_n = fxy ./ vGM_fx(i);
    gxy_n = gxy ./ vGM_gx(i);
    hxy_n = hxy ./ vGM_hx(i);
    
    % %
    % Preprocess
    [alpha, beta, th1, th2] = Preprocess_3Polys(fxy_n, gxy_n, hxy_n, m, n, o, k);
    
    vAlpha(i) = alpha;
    vBeta(i) = beta;
    vTh1(i) = th1;
    vTh2(i) = th2;
  
    % Get f(x,y) and g(x,y) with thetas to get f(w_{1},w_{2}) and g(w_{1},w_{2})
    fww = GetWithThetas(fxy_n, m, th1, th2);
    a_gww = alpha .* GetWithThetas(gxy_n, n, th1, th2);
    b_hww = beta .* GetWithThetas(hxy_n, o, th1, th2);
    
    
    % Build the Sylvester Matrix
    Sk = BuildSylvesterMatrix_3Polys(fww, a_gww, b_hww, m, n, o, k);
    
    % Get vector of singular values
    arr_SingularValues{i} = svd(Sk);
    
    [~,R] = qr(Sk);
    [~,c] = size(R);
    
    arr_R1{i} = R(1:c, 1:c);
    
end



switch SETTINGS.RANK_REVEALING_METRIC
    % R1 Row Norms
    % R1 Row Diagonals
    % Singular Values
    % Residuals
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
            
            plotRowNorms(arr_R1RowNorms, limits_k, limits_t)
            plotMaxMinRowNorms(vMaxRowNormR1, vMinRowNormR1, limits_k, limits_t, rank_range)
            
        end
        
        vMetric = log10(vMinRowNormR1 ./ vMaxRowNormR1);
        
    case 'R1 Row Diagonals'
        
        vMaxDiagonalEntry = zeros(nSubresultants, 1);
        vMinDiagonalEntry = zeros(nSubresultants, 1);
        
        for i = 1:1:length(arr_R1)
            
            % Get maximum diagonal
            vMaxDiagonalEntry(i) = max(abs(diag(arr_R1{i})));
            
            % Get minimum diagonal
            vMinDiagonalEntry(i) = min(abs(diag(arr_R1{i})));
            
        end
        
        vRatio_MaxMin_DiagonalEntry = vMinDiagonalEntry ./ vMaxDiagonalEntry;
        
        if(SETTINGS.PLOT_GRAPHS)
            
            %plotDiagonalsR1(arr_R1, myLimits, limits_t)
            plotMaxMinDiagonalR1(vRatio_MaxMin_DiagonalEntry, limits_k, limits_t, rank_range);
            
        end
        
        vMetric = log10(vRatio_MaxMin_DiagonalEntry);
        
    case 'Singular Values'
        
        vMinimumSingularValues = zeros(length(arr_SingularValues),1);
        for i = 1:1:length(arr_SingularValues)
            
            vec = arr_SingularValues{i};
            vMinimumSingularValues(i) = min(vec);
            
            
        end
        
        if(SETTINGS.PLOT_GRAPHS)
            
            plotSingularValues(arr_SingularValues, limits_k, limits_t);
            plotMinimumSingularValues(vMinimumSingularValues, limits_k, limits_t, rank_range);
            
        end
        
        vMetric = log10(vMinimumSingularValues);
        
    case 'Residuals'
        
        error('err')
        
    otherwise
        error('err')
end


if (lowerLimit_k == upperLimit_k)
    
    % Get the degree of the GCD
    t = GetGCDDegree_OneSubresultant(vMetric);
    
else
    
    % Get the degree of the GCD
    t = GetGCDDegree_MultipleSubresultants(vMetric, limits_k, rank_range);
    
    rank_range_low = vMetric(t - (lowerLimit_k - 1) );
    if (t < upperLimit_k)
        rank_range_high = vMetric(t - (lowerLimit_k - 1) + 1);
    end
    
    % update rank range
    rank_range = [rank_range_low rank_range_high];
   
    
    GM_fx = vGM_fx(t - lowerLimit_k + 1);
    GM_gx = vGM_gx(t - lowerLimit_k + 1);
    GM_hx = vGM_hx(t - lowerLimit_k + 1);
    alpha = vAlpha(t - lowerLimit_k + 1);
    th1 = vTh1(t - lowerLimit_k + 1);
    th2 = vTh2(t - lowerLimit_k + 1);
    
end

end