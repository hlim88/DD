function[kn] = stiffmatrix()
% The purpose of this function is to build the stiffness matrices
% for the Poisson equation on the 2-D unit square using triangular elements
% for both the Dirichlet and Neumann boundary conditions.
% The calling procedure is [kd, kn] = stiffmatrix(n) where
% n = 1/h, where he is the width of each element; 
% kd = (n-1)x(n-1) stiffness matrix for the free nodes under Dirichlet BC;
% kn = (n+1)x(n+1) stiffness matrix for the free nodes under Neumann BC. 
% This function only considers the case where h is equal in both the
% x and y directions.
global n;
global A;

y = n;
z = n + 2;

% Initialize matrices
A = zeros(y^2);
kn = zeros(z^2);

% Create Dirichlet matrix
% Diagonal block
w = diag(4*ones(y,1)) + diag(-1*ones(y-1,1),-1) + diag(-1*ones(y-1,1),1);

for i = 0:y-1
    A(i*y+1:(i+1)*y, i*y+1:(i+1)*y) = w;
    
    if i ~= 0
        A(i*y+1:(i+1)*y, (i-1)*y+1:i*y) = -eye(y);
    end
    
    if i ~= y-1
        A(i*y+1:(i+1)*y, (i+1)*y+1:(i+2)*y) = -eye(y);
    end
end

% Boundary diagonal block
db = diag(2*ones(z,1)) + diag(-0.5*ones(z-1,1),-1) + diag(-0.5*ones(z-1,1),1);
db(1, 1) = 1; db(z, z) = 1;

% Interior diagonal block
di = diag(2*ones(z,1)) + diag(-1*ones(z-1,1),-1) + diag(-1*ones(z-1,1),1);
di(2:z-1, 2:z-1) = A(1:y, 1:y);

% Off-diagonal block
od = diag(-0.5*ones(z,1));
od(2:z-1, 2:z-1) = -eye(y);

% Neumann stiffness matrix
for i = 0:z-1
    if i == 0 || i == z-1
        kn(i*z+1:(i+1)*z, i*z+1:(i+1)*z) = db;
    
    else
        kn(i*z+1:(i+1)*z, i*z+1:(i+1)*z) = di;
    end
    
    if i ~= 0
        kn(i*z+1:(i+1)*z, (i-1)*z+1:i*z) = od;
    end
    
    if i ~= z-1
        kn(i*z+1:(i+1)*z, (i+1)*z+1:(i+2)*z) = od;
    end
end

return