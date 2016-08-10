function arr_hx = Deconvolve_Bivariate_Separate(arr_fx)
% Peform the sequence of deconvolutions independently

% Get number of polynomials f_{i}(x)
nPolys = length(arr_fx);

% Initialise an array to store polynomials h_{i}(x) = h{i-1} ./ h{i}
arr_hx = cell(nPolys-1,1);

% For each pair of polynomials perform deconvolution.
for i = 1:1:nPolys-1
    
    arr_hx{i} = Deconvolve_Bivariate(arr_fx{i},arr_fx{i+1}) ;
    
end

end