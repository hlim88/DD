global n;   n = 7;
global nc;  nc = 3;

for i = 1:n*n
    x(i, 1) = i;
end

[R] = restrict(n, nc);

%y = R * x;

%z = R' * y;

%x
%y
%z


