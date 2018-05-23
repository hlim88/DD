n = 7;
ovlp = 2;
[A, An] = stiffmatrix(n+1);
[ind_a11, ind_a12, ind_a21, ind_a22] = index_sub(n, ovlp);

m = (n+1)/2;
for i = 1:(m+ovlp)^2
    for j = 1:(m+ovlp)^2
        Asub(i,j) = A(ind_a11(i), ind_a11(j));
    end
end

[L, U, P, flag] = LUGEN(Asub);

y = zeros(n^2,1);
