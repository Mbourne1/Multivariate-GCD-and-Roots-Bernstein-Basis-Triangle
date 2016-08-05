function Q = BuildQ(m,n,k)

Q1 = BuildQ1(n-k);
Q2 = BuildQ1(m-k);

Q = blkdiag(Q1,Q2);

end