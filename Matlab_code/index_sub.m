function index_sub()

% The purpose of this function is to create vectors that relate the
% local node numbers on the subdomains to the global node numbers.
% The calling procedure is [w, x, y, z] = index_sub(n, overlap) where:
% n = number of nodes in a row of the global domain;
% overlap = number of nodes past the middle node of the global domain;
% w = vector of global node numbers for subdomain A11;
% x = vector of global node numbers for subdomain A12;
% y = vector of global node numbers for subdomain A21; and
% z = vector of global node numbers for subdomain A22.

global n; global ind_a11; global ind_a12; global ind_a21; global ind_a22;
global ovlp;

% Determine middle node of the global domain
m = (n+1)/2;

% Assign global node numbers in the local vectors
k = 1;
for i = 1:m+ovlp
    for j = 1:m+ovlp
        ind_a11(k) = j + (i-1)*n;
        ind_a12(k) = m - ovlp + (j - 1) + (i - 1)*n;
        ind_a21(k) = (m - ovlp + (i - 2))*n + j;
        ind_a22(k) = (m - ovlp + (i - 2))*n + (j - 1) + m - ovlp;
        k = k + 1;
    end
end

return