function [fxy_noisy,noise] = Noise2(fxy,el)

fxy_noisy = fxy;
[m1,m2] = GetDegree(fxy);

noise = zeros(m1+1,m2+1);


end