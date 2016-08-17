function Cf = BuildT1(fxy,m,n_t)
% Build the matrix C where C(f(x,y))*v = h. 

% Inputs.
%
% fxy : coefficients of polynomial f(x,y)
%
% m : Total degree of f(x,y)
%
% n-t : Total degree of v(x,y)

% Get number of coefficients of nchoosek
nCoefficients_gxy = nchoosek(n_t+2,2);

% Get number of coefficients in the product fg = h(x,y)
nCoefficients_hxy = nchoosek(m+n_t+2,2);

% Initialise a zero matrix
zero_mat = zeros(m+n_t+1,m+n_t+1);

Cf = zeros(nCoefficients_hxy,nCoefficients_gxy);

% Get fxy with trinomial coefficients
fxy_tri = GetWithTrinomials(fxy,m);


count = 1;

for diag_index = 0:1:n_t
    for i = diag_index:-1:0
        j = diag_index - i;
        
        % Get matrix of coefficients of f(x,y) shifted down by i rows and
        % across by j columns
        temp_mat = zero_mat;
        temp_mat(i+1:i+m+1,j+1:j+m+1) = fxy_tri;
        
        temp_vec = GetAsVector(temp_mat);
        
        % Remove all but the first nchoosek(m+n+2,2) coefficients
        temp_vec = temp_vec(1:nCoefficients_hxy);
        
        % Insert coefficients into the i,j th column of C(f(x,y)).
        
        Cf(:,count) = temp_vec;
        
        % Increment counter
        count = count + 1;
        
    end
end



end