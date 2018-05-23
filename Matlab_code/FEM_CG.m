function main
% The purpose of this function is to solve the system Ax = b
% for an FEM formulated Dirichlet stiffness matrix using CG method.

clear all

% User chosen values
max_it = 100;
global n; n = 51;              %number of nodes in a row of Dirichlet
global ovlp; ovlp = 2;
tol = 1e-6;

% Set up global variables
global A; global L; global U; global P;
global ind_a11; global ind_a12; global ind_a21; global ind_a22;

% Formulate a solution: sin(pi*x)*sin(pi*y)
for i = 1:n
    for j = 1:n
        s((i-1)*n + j) = sin(pi*(i-1)/(n-1))*sin(pi*(j-1)/(n-1));
    end
end

% System matrix
[An] = stiffmatrix();

% Right-hand side
b = A * s';

% Initial guess vector
x = zeros(n^2, 1);

% Define index vectors
index_sub();

% Form domain submatrices
m = (n+1)/2;
for i = 1:(m+ovlp)^2
    for j = 1:(m+ovlp)^2
        Asub(i,j) = A(ind_a11(i), ind_a11(j));
    end
end

% Compute PA = LU factorization on the submatrices
[flag] = LUGEN(Asub);

% CG solver
[y, error, iter, flag] = cgm(x, b, max_it, tol);

if flag == 1
    disp('Method did not converge')
end

% Output plot
for i = 1:n
    for j = 1:n
        y_surf(i, j) = y(i + (j-1)*n);
        s_surf(i, j) = s(i + (j-1)*n);
    end
end

subplot(2,2,1); surf(y_surf)
title('Estimated')
zlim([0,1])
subplot(2,2,2); surf(s_surf)
title('Actual')
zlim([0,1])
subplot(2,2,3); surf(real(s_surf - y_surf))
title('Difference')



