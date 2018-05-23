function [z] = addsm(r)

global n; global U; global L; global P;
global ind_a11; global ind_a12; global ind_a21; global ind_a22;

z = zeros(n^2, 1);


for i = 1:2
    for j = 1:2
        q = ['z',int2str(i),int2str(j),'= U \(L \ P) * r(ind_a',int2str(i), int2str(j), ')'];
        eval(q);
       
        w = ['z(ind_a', int2str(i), int2str(j), ') = z(ind_a', int2str(i), int2str(j), ') + z',int2str(i),int2str(j)];
        eval(w);
   end    
end

return