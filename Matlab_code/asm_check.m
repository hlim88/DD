function main
n = 5;      %nodes in a row
ovlp = 1;
m = (n + 1)/2;
U = eye((m+ovlp)^2);
L = eye((m+ovlp)^2);
P = eye((m+ovlp)^2);

for i = 1:n^2
    g(i) = i;
end
r = g';
[ind_a11, ind_a12, ind_a21, ind_a22] = index_sub(n, ovlp);

z = addsm(n, U, L, P, r, ind_a11, ind_a12, ind_a21, ind_a22)

return