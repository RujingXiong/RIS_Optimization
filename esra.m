function [w] = esra(A,B)
%求解max(||Aw||2),B为bit数
[~,N] = size(A);
omega = 2*pi/(2^B);
w = exp(1i*omega.*randi(2^B,N,1));
deta = [0,1];
deta(1) = norm(A*w);
while(deta(2) > 0)
    for i = 1:N
        w1 = w;
        w1(i,1) = 0;
        aaa = w1'*(A')*A*w1;
        w1(i,1) = 1;
        bbb = w1'*(A')*A*w1 - aaa;
        for j = 1:2^B
            if distance(angle(bbb)+j*omega,angle(aaa)) <= omega/2
                w(i,1) = exp(1i*omega*j);
            end
        end
    end
    deta(2) = abs(deta(1)-norm(A*w));
    deta(1) = norm(A*w);
end

