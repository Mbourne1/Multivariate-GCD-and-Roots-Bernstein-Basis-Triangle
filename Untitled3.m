m = 7;
n = 12;
k = 4;


exp1 = 3 ./ (nchoosek(m + 2, 2) * nchoosek(n - k + 1, 2));
exp2 = 3 ./ (nchoosek(m + 2, 2) * nchoosek(n - k + 2, 2));

expression_a = exp1 - exp2;


%
%
%
%

expression_b = 2 / (nchoosek(m + 2,2) * nchoosek(n - k + 2,3));


display(expression_a)
display(expression_b)