%%%%--------------------Compare with Exhaustive---------------------
%Author£ºJialong Lu, Rujing Xiong
%Email: jialong@hust.edu.cn, rujing@hust.edu.ch
% Revise History:
%  Created 
% 2024.05.27  Modified 
% Without los
% Description£º Exhaustive_search
clear all
clc
tic
T = 3;
N = 5;
sigma2 = 1;
B=1; % number of bits
y_ = zeros(1,length(N));
y_2 = zeros(1,length(N));
loop = 1000;
y_square = zeros(loop,length(N));
y_square2 = zeros(loop,length(N));
for m = 1:loop    
    for t = 1:length(N)    
        G = (randn(N(t),T)+1i*randn(N(t),T)).*sqrt(1/2);     %AP-RIS  
        Hr = (randn(1,N(t))+1i*randn(1,N(t))).*sqrt(1/2);    %RIS-USER
      % Hr = [3*ones(N(t)/3,1); 1*ones(N(t)/3,1);5*ones(N(t)/3,1)];%RIS-USER      
        Phi = diag(Hr)*G;
      % R = Phi*Phi';
       %% Exhaustive search   
        w = 0:1:2^N(t)-1;
        w = dec2bin(w,N(t));
        W = (w == '1');
        W = double(W);
        W(W==0) = -1;
        Y = W*Phi;
        norm_Y = vecnorm(Y,2,2); %2-norm of each row in Y

        %norm_Y = diag(norm_Y);% the elements on the diagonal.
        [y_(t),w_index] = max(norm_Y);
        w_opt_e = W(w_index,:);
        
        
        % Discret optimization-lowrank
        Vector1 = finds2(Phi,B); % N>=2T-1
        Y2 = Vector1.'*Phi;
        norm_Y2 = vecnorm(Y2,2,2); %2-norm of each row in Y2
        %Y = W'*R*W;
        [y_2(t),w_index2] = max(norm_Y2);
        w_opt = Vector1(:,w_index2);


    end
    y_square(m,:) =y_.^2;
    y_square2(m,:) = y_2.^2;
    X = sprintf('The loop have completed %d times.',m);
    disp(X);
end
y_e = mean(y_square,1);           
SNR_e = 10*log10(y_e.^2./sigma2); 
y = mean(y_square2,1);          %
SNR = 10*log10(y.^2./sigma2); 
toc
figure
plot(N,SNR,'-rs','LineWidth',2);
hold on
plot(N,SNR_e,'b*','LineWidth',3);
xlabel('N','Interpreter','latex','Fontsize',15);
ylabel('SNR(dB)','Interpreter','latex','Fontsize',15);
legend('Proposed Method','Exhaustive Search')
grid on
title('');
toc