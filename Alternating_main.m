%% Author:  Rujing Xiong, Tiebin Mi
 % CreatTime: 2024.02.20
 % Complete:
 % Modified: 
 % E-mail: Rujing@hust.edu.cn
 % Description: Alternating inner product maximization approach, mainly for
 % p-norm ||Aw||_p(p is 2 within this code), while A is a matrix with
 % rank-M, but not rank-1
 % 
 % 
clear all
clc
tic

N = 100;
M = 16;
numBits = 1;

A = randn(M,N)+1i*randn(M,N);

epli = 1;
R = A'*A;

%%%%%%%%%%%%%%  Proposed Alternating Maximization   %%%%%%%%%%%%%%
W = exp(1i*randn(N,1)); % W can be selected randomly, also can be selected through some optimization methods, such as a quantized version from continuous solution.
while epli>1e-10
    fval_1 = norm(A*W,2);
    Z = ((A*W)/norm(A*W,2));  %  Z为M*1
    W = (A'*Z)/norm(A'*Z);
    fval_2 = norm (A*W,2);
    epli = abs(fval_2-fval_1);
end
disp(fval_2);
epli = 1; % recover epli for the proceeding of next loop

while epli>1e-10
    fval_1 = norm(A*W,2);
    Z = ((A*W)/norm(A*W,2));  %  Z为M*1
    Temp = A'*Z;
    W = Opt_absWZ(Temp,numBits);
    fval_2 = norm (A*W,2);
    epli = abs(fval_2-fval_1);
end
disp(fval_2);
y_my = fval_2;
w_my = W;

%%%%%%%%%%%%%%%   Random      %%%%%%%%%%%%%%%%%%%%%%%%%%
rand = 100000;
w_rand = exp(1i*randn(N,rand));
Y = A*w_rand;
column_norms = vecnorm(Y);
[y_rand,index] = max(column_norms);
% ww = w_rand(:,index);
% ww = Discretization2Bit(ww);
% y_value2 =  norm(A*ww);
disp(y_rand);

%%%%%%%%%%%%%%%%    Manopt    %%%%%%%%%%%%%%%%%%%%%%%%%%     
manifold = complexcirclefactory(N);
problem.M = manifold;
problem.cost = @(w) -w'*R*w; %
problem.grad = @(w) manifold.egrad2rgrad(w,-2*R*w);
[w_man,wcost,info,options] = steepestdescent(problem);
%w_man = Discretization2Bit(w_man);
w_man = discrete(w_man,numBits,N);
y_man = norm(A*w_man);

%%%%%%%%%%%%%%%%   SDR  %%%%%%%%%%%%%%%%%%%%%%%%%%%
w_sdr = sdr(N-1,R,1);
w_sdr = discrete(w_sdr,numBits,N);
y_sdr = norm(A*w_sdr);

%%%%%%%%%%%%%%%%   Successive refinement  %%%%%%%%%
        w_esra = esra(A,numBits);
        y_esra = norm(A*w_esra);
%%%%%%%%%%%%%%%% Exhaustive search   %%%%%%%%%%%%%%%%%%%%
w_ex = 0:1:2^N-1;
w_ex = dec2bin(w_ex,N);
W_ex = (w_ex == '1');
W_ex = double(W_ex);
W_ex(W_ex==0) = -1;
Y_ex = A*(W_ex.');
norm_Y = vecnorm(Y_ex); %2-norm of each row in Y
%norm_Y = diag(norm_Y);%取对角线上的每一个元素
[y_ex,w_index] = max(norm_Y);
w_ex = W_ex(w_index,:);
toc
