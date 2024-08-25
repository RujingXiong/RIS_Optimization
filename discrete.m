function [opt] = discrete(w,B,n)
for i=1:1:n
    for j=0:1:2^B-1
        if distance(angle(w(i,1)),j*2*pi/2^B) < pi/2^B
           w(i,1) = exp(1i*j*pi/2^(B-1));
        end
    end
end
opt = w;