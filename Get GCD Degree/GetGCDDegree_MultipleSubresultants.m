function [t] = GetGCDDegree_MultipleSubresultants(vMetric, limits_k, rank_range)
%
% % Inputs
%
% vMetric : (Vector) 
%
% limits_k : (Int Int)
%
% rank_range : [Float Float]

global SETTINGS

% Get the name of the function which called this function
[St,~] = dbstack();
calling_function = St(2).name;

% Get maximum change in the rank revealing metric
[maxDelta, index] = Analysis(vMetric);

% Get upper and lower bounds
lowerLimit_k = limits_k(1);
upperLimit_k = limits_k(2);

% check if the maximum change is significant
fprintf([mfilename ' : ' sprintf('Max change : %2.4f \n', maxDelta)]);
fprintf([mfilename ' : ' sprintf('Previous Max Change : %2.4f \n', abs(diff(rank_range)))]);



if abs(maxDelta) < (0.55 * abs(diff(rank_range)))
    fprintf([calling_function ' : ' mfilename ' : ' 'Polynomials either coprime or GCD = g(x) \n' ])
    
    
    % Change in Singular values is not significant so check if all
    % subresultants are rank deficient or full rank
    
    % Get the average of the minimum singular values
    avg = mean(vMetric);
    
    fprintf([mfilename ' : ' sprintf('THRESHOLD RANK :  %2.4e \n', SETTINGS.THRESHOLD_RANK)]);
    fprintf([calling_function ' : ' mfilename ' : ' sprintf('Average Singular Value : %e \n',avg) ])
    
    if avg > mean(rank_range)
       % All Minimum singular values are below threshold so, all 
       % subresultants are rank deficient. deg(GCD) = 0
       fprintf([calling_function ' : ' mfilename ' : ' 'Polynomails are coprime\n' ])
       t = 0;
       
    else 
        % All minimum singular values are above threshold so all
        % subresultants are full rank. deg(GCD) = min(m,n)
       fprintf([calling_function ' : ' mfilename ' : ' 'All Subresultants are rank deficient : GCD = g(x) \n' ])
       t = upperLimit_k;
       
    end

else
    % change is significant
    t = lowerLimit_k + index - 1;

    
end

end