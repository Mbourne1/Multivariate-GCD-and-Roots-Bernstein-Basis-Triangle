function DT1Q1 = BuildDT1Q1(fxy, m, n_k)


D = BuildD_2Polys(m,n_k);

Q1 = BuildQ1(n_k);

T1 = BuildT1(fxy,m,n_k);

DT1Q1 = D*T1*Q1;

end