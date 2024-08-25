function w_sdr = sdr(N,R,count)
f_tmp = 0;
r_tmp = zeros(N,1);   % r=[];
w_tmp = zeros(N,1);
for k=1:count
    r = (randn(N+1,1)+1i*randn(N+1,1)).*sqrt(1/2);   % (N,1)
    cvx_begin
    variable V(N+1,N+1) symmetric semidefinite   %变量是一个(N)*(N)的对称半正定矩阵
    maximize( real(trace(R*V)))
    subject to
    diag(V) == 1;
    cvx_end

    [U,Sigma] = eig(V);
    w = U*Sigma^(1/2)*r;   % (N*1)
    f = w'*R*w;       %随机次数为count次，找到其中最大的f对应的波束赋形向量w和高斯随机向量r
    if f>f_tmp
        f_tmp = max(f,f_tmp);
        r_tmp = r;
        w_tmp = w;     %求解出来的w_tmp为啥比w_的维度要小10的倍数个元素
    end
end
%     [m,index]=max(f);
%     r_opt = r(:,index);
%     w_opt = w(:,index);
%     theta_opt = angle(w_opt);
%     w_opt = exp(1i*theta_opt);     % 使其满足恒模约束
%     W = diag(w_opt);
w_tmp = w_tmp./w_tmp(N+1);
theta_opt = angle(w_tmp);
w_sdr = exp(1i*theta_opt);