function [R] = restrict (n, nc);
%global n;
%global nc;

R = zeros(nc*nc, n*n);
hc = 1 / (nc - 1);
h = 1 / (n - 1);

for i = 1:nc
    yc = (i - 1) * hc;
    
    for j = 1:nc
        xc = (j - 1) * hc;
        
        for k = 1:n
            y = (k - 1) * h;
                
            for m = 1:n
                x = (m - 1) * h;
                [w] = basis(x, y, xc, yc, hc);
                R(j+(i-1)*nc, m+(k-1)*n) = w;
            end
        end
    end
end

