function w = huber(x, r)
[m, n] = size(x);

for i = 1:m
    for j = 1:n
        if abs(x(i,j))<r
            w(i,j) = x(i,j)^2/2;
        else
            w(i,j) = r*abs(x(i,j)) - r^2/2;
        end                
    end
end