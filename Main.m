%% Author: Rujing Xiong, Tiebin Mi
 % CreatTime: 2022.10.07
 % Complete:2022.10.07
 % Modified:
 % Modified:
 % E-mail: Rujing@hust.edu.cn
 % Description: DiscreteOptimization for 2-bit RIS, compare with Manopt, for multi-bit
 % Based on traditional channel model, Mi's code
 %  SNR with N(t) 2bit and multi bit
 clc
 clear
 tic
 sigma2 = 1;  %noise
 T = 1;  %number of Tx
 B = 2; % number of bits
 N = 100:100:1000;

 y_square = zeros(1,length(N));

 loop = 1000;
 y_my = zeros(loop,length(N)); % my method (proposed DaS)
 y_man = zeros(loop,length(N)); % Manopt-discrete
 y_Man = zeros(loop,length(N)); % Manopt
 y_apk = zeros(loop,length(N)); % apk
 y_sdr = zeros(loop,length(N)); % SDR-SDP-discrete

 for h = 1:loop
     for t = 1:length(N)
         U = zeros(N(t),2^B*N(t));
         H = (randn(N(t),T)+1i*randn(N(t),1)).*0.5;     %AP-RIS
         Hr = (randn(N(t),1)+1i*randn(N(t),1));%RIS-USER
         Phi = diag(Hr)*H;
         R = Phi*Phi';

         [Z,eigs_R] = eig(R,'vector');
         [eigs_R,index_z] = sort(eigs_R,'descend');
         Z = Z(:,index_z);  %The corresponding eigenvectors are sorted in descending order too.
         Z_S = Z(:,1);      %Taking the igenvector corresponding to the largest eigenvalue
         
         % proposed DaS
         [w_opt,opt_my] = Opt_absWZ(Z_S,B);
         y_my(h,t) = w_opt'*R*w_opt;
        
         %  Manopt
         manifold = complexcirclefactory(N(t));
         problem.M = manifold;
         problem.cost = @(w) -w'*R*w; %
         problem.grad = @(w) manifold.egrad2rgrad(w,-2*R*w);
         [w,wcost,info,options] = steepestdescent(problem); % at a random point on the manifold
         %[w,wcost,info,options] = conjugategradient(problem);%
         %[w,wcost,info,options] = barzilaiborwein(problem);  %Barzilai Borwein
         %y_Man(h,t) = w'*R*w; % under continue phase shift
         %w = Discretization2Bit(w);
         w = discretization1Bit(w);
         y_man(h,t) =w'*R*w; % under discrete phase shift

         % APX
         [w_apk,opt_apk] = APK_K_Ary(Z_S,B);
         y_apk(h,t) = w_apk'*R*w_apk;
         
         % CVX  
         f_tmp = 0;
         r_tmp = zeros(N,1);   % r=[];这三个tmp的初始化不能放在循环外面，应该要在每当N变化的时候就要初始化，这样的话才能使得后面对于r_tmp.w_tmp的赋值不会出现延时。
         w_tmp = zeros(N,1);
         count = 10;
         for k=1:count
             r = (randn(N,1)+1i*randn(N,1)).*sqrt(1/2);   % (N,1)
             cvx_begin
             variable V(N,N) symmetric semidefinite   %变量是一个(N)*(N)的对称半正定矩阵
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
         theta_opt = angle(w_tmp);
         w_sdr = exp(1i*theta_opt);
         w_sdr = discretization1Bit(w_sdr);
         y_sdr(h,t) = w_sdr'*R*w_sdr;
     end
     X = sprintf('The loop have completed %d times.',h);
     disp(X);
 end % end loop

 y_meanMy = mean(y_my,1);
 SNR_my = 10*log10(1+y_meanMy./sigma2);

 y_meanMan = mean(y_Man,1);
 SNR_Man = 10*log10(1+y_meanMan./sigma2);

 y_meanman = mean(y_man,1);
 SNR_man = 10*log10(1+y_meanman./sigma2);

 y_meanapk = mean(y_apk,1);
 SNR_apk = 10*log10(1+y_meanapk./sigma2);

 y_meansdr = mean(y_sdr,1);
 SNR_sdr = 10*log10(1+y_meansdr./sigma2);

figure
 hold on
 plot(N,SNR_Man,'-rs','LineWidth',1);
 hold on
 plot(N,SNR_man,'-gs','LineWidth',1);
 hold on
 plot(N,SNR_my,'--cd','LineWidth',1);
 hold on
 plot(N,SNR_apk,'--bx','LineWidth',1);
 hold on 
 plot(N,SNR_sdr,'--gd','LineWidth',1);
 xlabel('N','Interpreter','latex','Fontsize',15);
 ylabel('SNR(dB)','Interpreter','latex','Fontsize',15);
 legend('Manopt','Discrete Manopt','Proposed method','APK','SDR-SDP','Interpreter','latex','Fontsize',15);
 hold on ;
 grid on;
 hold on;
 box on;
% title('Gains on Signal power');
 toc

% figure
%  hold on
%  plot(N,SNR1,'-gs','LineWidth',1);
%  hold on
%  plot(N,SNR2,'--cd','LineWidth',1);
%  hold on
%  plot(N,SNR3,'-bx','LineWidth',1);
%   hold on
%  plot(N,SNR4,'-bs','LineWidth',1);
%   hold on
%  plot(N,SNR5,'--gd','LineWidth',1);
%   hold on
%  plot(N,SNR6,'-cx','LineWidth',1);
%   hold on
%  plot(N,SNR7,'-g+','LineWidth',1);
%   hold on
%  plot(N,SNR8,'--b+','LineWidth',1);
%  hold on
%  plot(N,SNR0,'-rs','LineWidth',1);
%  xlabel('N','Interpreter','latex','Fontsize',15);
%  ylabel('SNR(dB)','Interpreter','latex','Fontsize',15);
%  legend('1-bit','2-bit','3-bit','4-bit','5-bit','6-bit','7-bit','8-bit','Manopt','Interpreter','latex','Fontsize',15);
%  hold on ;
%  grid on;
%  hold on;
%  box on;
%  toc

% figure
%  hold on
%  plot(N,SNR1,'-gs','LineWidth',1);
%  hold on
%  plot(N,SNR2,'--cd','LineWidth',1);
%  hold on
%  plot(N,SNR3,'-bx','LineWidth',1);
%   hold on
%  plot(N,SNR4,'-bs','LineWidth',1);
%   hold on
%  plot(N,SNR5,'--gd','LineWidth',1);
%  hold on
%  plot(N,SNR0,'-rs','LineWidth',1);
%  xlabel('N','Interpreter','latex','Fontsize',15);
%  ylabel('SNR(dB)','Interpreter','latex','Fontsize',15);
%  legend('1-bit','2-bit','3-bit','4-bit','5-bit','Manopt','Interpreter','latex','Fontsize',15);
%  hold on ;
%  grid on;
%  hold on;
%  box on;
%  toc

%set(gca,'yticklabel',{'5%','6%','7%','8%','9%','10%'})%set ylabel in 100%
