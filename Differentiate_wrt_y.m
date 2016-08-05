function [d_fxy] = Differentiate_wrt_y(fxy)



fxy = transpose(fxy);

d_fxy = Differentiate_wrt_x(fxy);

d_fxy = transpose(d_fxy);


end