function [flag] = LUGEN(w)
% The purpose of this function is to explicitly find the LU factoriztion of 
% a square matrix using Gaussian Elimination with partial pivoting.
% The calling procedure is [L, U, p, flag] = LUGEN(w) where:
% w = square coefficient matrix; L = lower triangular matrix;
% U = upper triangular matrix; P = permutation matrix;
% flag = 0 if a solution is found or 1 if there is no unique solution.

global L; global U; global P;

[m,n] = size(w);
flag = 0;

fac = w;
L = eye(n);
P = zeros(n);
U = zeros(n);

if m~=n
    disp('The input matrix must be a square matrix')
    return
end

for i=1:n
    q(i) = i;
end

for i=1:n-1
    
    %Find the pivot element and check to make sure it is nonzero.
    
    [facmax, pivot]=max(abs(fac(i:n,i)));
    if facmax == 0
        flag = 1;
        return
    end
    
    pivot = pivot + i - 1;
    
    %Interchange rows if necessary
    
    if pivot ~= i
        temp = fac(i,1:n);
        fac(i,1:n) = fac(pivot,1:n);
        fac(pivot,1:n) = temp;
        
        %Interchange rows of vector of row interchanges
        
        temp = q(i);
        q(i)= q(pivot);
        q(pivot) = temp;
    end
    
    %Apply the elementary row operations to zero the subdiagonal entries
    %in column i, while storing the multiplier where the zero would go.
    
    for j=i+1:n
        m = -fac(j,i) / fac(i,i);
        fac(j,i) = -m;
        fac(j,i+1:n) = fac(j,i+1:n) + m * fac(i,i+1:n);
    end
end

if fac(n,n) == 0
    flag = n;
    return
end

for i = 1:n
    for j = 1:n
        if i <= j
            U(i,j) = fac(i,j);
        else
            L(i,j) = fac(i,j);
        end
    end

    P(i, q(i)) = 1;
end