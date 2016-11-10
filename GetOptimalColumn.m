function idx_col = GetOptimalColumn(Sk)
%% Find Optimal column for removal from S_{t_{1},t_{2}}
% Given that t1 and t2 have been calculated build the sylvester matrix and
% find the optimal column such that a residual is minimized
%
% % Inputs
%
% Sk : The Sylvester subresultant matrix S_{k}(f,g)
%
% % Outputs
%
% idx_col : Index of optimal column for removal

% Get the number
[~,nCols_Sk] = size(Sk);

% QR Decomposition of the Sylvester Matrix S_{k}
[Qk,Rk] = qr(Sk);

% For each column of S_{t_{1},t_{2}}, move the column to the RHS, to obtain
% Ax=b, where A consists of the remaining columns of S_{t_{1},t_{2}} and b
% is the removed column c_{k}

% Intialise a zero vector to store residuals obtained by the removal of
% each column.
vec_residuals_QR = zeros(nCols_Sk,1);

for k = 1 : 1 : nCols_Sk
    
    % Get column c_{k} for removal from S_{k}
    ck = Sk(:,k);
    
    % Perform QR delete to remove k column from QR decomposition of 
    % S_{t_{1},t_{2}}
    [Q,~] = qrdelete(Qk,Rk,k);
    
    cd = Q'*ck;
    
    d = cd(nCols_Sk+1:end,:);
    
    % Get Residuals
    vec_residuals_QR(k) = norm(d);
    
end

%Obtain the column for which the residual is minimal.
[~,idx_col] = min(log10(vec_residuals_QR));

% Print out optimal column for removal.
fprintf([mfilename ' : ' sprintf('Optimal column for removal is : %i \n',idx_col)]);

end