syms x y


f0 = (x - 1)^3*(x - 2)^5*(y - 3)^5;

m = double(feval(symengine, 'degree', f0));
fxy = double(rot90(coeffs(f0,[x,y],'All'),2));

fxy_big = zeros(m+1,m+1);
[nRows,nCols] = size(fxy);
fxy_big(1:nRows,1:nCols) = fxy;

fxy_bb = PowerToBernstein(fxy_big,m);



% Get derivative of f(x,y) in power basis
df_dx = diff(f0,x)

dfxy_dx = double(rot90(coeffs(df_dx,[x,y],'All'),2));
dfxy_dx_big = zeros(m,m)
[nRows,nCols] = size(dfxy_dx)
dfxy_dx_big(1:nRows,1:nCols) = dfxy_dx;



% Get the derivative of f(x,y) in the Bernstein basis
dfxy_bb_dx = Differentiate_wrt_x(fxy_bb)

% Check that derivative in power basis converts to same derivative in
% bernstein basis
dfxy_dx_converted = PowerToBernstein(dfxy_dx_big,m-1)
dfxy_bb_dx./dfxy_dx_converted