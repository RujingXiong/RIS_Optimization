%%%%--------------------Exhaustive Search ---------------------
%Author£ºRujing Xiong
%Email:Rujing@hust.edu.ch
% Revise History:
% 2022-06-09  Created 
% 
% with los
% Description£º A comparison of the signal-to-noise ratio (SNR) performance between the exhaustive search and the proposed DaS algorithm
clear all
clc
tic
T = 1;
N = 1:1:10;
sigma2 = 1;

y_square_e = zeros(1,length(N));
y_square = zeros(1,length(N));
loop = 10000;
y__e = zeros(loop,length(N));
y_ = zeros(loop,length(N));
for m = 1:loop    
    for t = 1:length(N) 
        H = (randn(N(t),T)+1i*randn(N(t),T)).*2;     %AP-RIS  
        Hr = (randn(N(t),1)+1i*randn(N(t),1));%RIS-USER
        G = (randn(1,T)+1i*randn(1,T)).*sqrt(1/2);  %LOS
       
        Phi = diag(Hr)*H;
        %R = Phi*Phi';
        R = [Phi*Phi',Phi*G';G*Phi',G*G'];
        
     %% Exhaustive search   
        w = 0:1:2^N(t)-1;
        w = dec2bin(w,N(t));  %decimalism to bin
        W = (w == '1');     %logic value
        W = double(W);
        W(W==0) = -1;

        W_last = ones(2^(N(t)),1);  % The last element replenished of each row in W
        W_bar = [W,W_last];
        
        
        Y = W_bar*R*W_bar';
        Y = diag(Y);
        [y_square_e(t),w_index] = max(Y);
        w_opt_e = W(w_index,:);
        w_opt_efinal = w_opt_e(1,N(t));

        %% Discret optimization(DaS)
        [Z,eigs_R] = eig(R,'vector');
        [eigs_R,index_z] = sort(eigs_R,'descend');
        Z = Z(:,index_z);  %The corresponding eigenvectors are sorted in descending order too.
        Z_S = Z(:,1);      %Taking the eigenvector corresponding to the largest eigenvalue
        %
        Phi_i = angle(Z_S);
        Phi_i(Phi_i>-pi & Phi_i<-pi/2) = Phi_i(Phi_i>-pi & Phi_i<-pi/2)+2*pi;
        I_1 = find(Phi_i>-pi/2 & Phi_i<pi/2);
        I_2 = find(Phi_i>pi/2 & Phi_i<3*pi/2);
        Phi_i_hat = Phi_i;
        Phi_i_hat(I_2) = Phi_i_hat(I_2)-pi;
        [Phi_i_hat,i_index] = sort(Phi_i_hat); % sort by increasing order
        %
        %
        W_hat = ones(N(t)+1,N(t)+1);
        for k = 1:N(t)+1
            W_hat(k,i_index(k+1:N(t)+1)) = -1;
        end
        for j = 1:N(t)+1
            W_hat(j,I_2) = -1*W_hat(j,I_2);
        end
        W = W_hat;
        W = W.';% each w in W is a column vector
        for e = 1:N(t)+1
            W(:,e) = W(:,e)./W(N(t)+1,e);%each w keeps the last element 1
        end
        Y = W'*R*W;
        Y = diag(Y);% get the diagonal elements
        [y_square(t),w_index] = max(Y);
        w_opt = W(:,w_index);
        w_my = w_opt(1:N(t));
    end
    y__e(m,:) = y_square_e;
    y_(m,:) = y_square;
    X = sprintf('The loop have completed %d times.',m);
    disp(X);
end

y_e = mean(y__e,1);          
SNR_e = 10*log10(y_e./sigma2); 

y = mean(y_,1);          %
SNR = 10*log10(y./sigma2); 
toc
figure
plot(N,SNR,'-rs','LineWidth',2);
hold on
plot(N,SNR_e,'b*','LineWidth',3);
xlabel('N','Interpreter','latex','Fontsize',15);
ylabel('SNR(dB)','Interpreter','latex','Fontsize',15);
legend('Proposed DaS','Eexhaustive Search')
grid on
title('');
toc