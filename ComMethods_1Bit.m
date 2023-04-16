%% author：Rujing Xiong
% time： 2022.6.05 created
% description ：Comparison of SNR Performance among Different Methods, including APX, Man, SDR, DaS 
% with los
%linear array

clear all
clc; tic

sigma2 = 1;  % noise。
T = 1; % Transmitter antennas
N = 10:5:50;
y_square_APX = zeros(1,length(N));
y_square_Man = zeros(1,length(N));
y_square_SDR = zeros(1,length(N));
y_square_DaS = zeros(1,length(N));


loop = 1000;
count = 1;
y_APX = zeros(loop,length(N));
y_Man = zeros(loop,length(N));
y_SDR = zeros(loop,length(N));
y_DaS = zeros(loop,length(N));

for m = 1:loop
    for t = 1:length(N)
        H = (randn(N(t),T)+1i*randn(N(t),T)).*2;     %AP-RIS
        Hr = (randn(N(t),1)+1i*randn(N(t),1));%RIS-USER
        G = (randn(1,T)+1i*randn(1,T)).*sqrt(1/2);  %LOS
%       G = zeros(1,T);

        Phi = diag(Hr)*H;

        R = [Phi*Phi',Phi*G';G*Phi',G*G'];

        %%% APX
        w_m = zeros(2,2*N(t)+2);
        y_star= zeros(2,2*N(t)+2);
        x_star = zeros(1,N(t));
        g_ystar = zeros(1,2*N(t)+2);
        V_V = zeros(1,N(t)+1);

        h_m = [diag(Hr)*H;G]';
        Re_hm = real(h_m);
        Im_hm = imag(h_m);
        V = [Re_hm;Im_hm];
        v_0 = V(:,1);

        R1 = h_m'*h_m;

        Phi_ = angle(h_m);
        y_minus =  exp(1i*(Phi_-pi/2));
        y_plus = exp(1i*(Phi_+pi/2));

        angle_ym = angle(y_minus);
        angle_yp = angle(y_plus);
        angle_all = [angle_ym,angle_yp];

        y_pie = exp(1i*angle_all+0.001); % all y' which will move in the interior
        re_ypie = real(y_pie);
        im_ypie = imag(y_pie);
        y_bold = [re_ypie;im_ypie];

        y_endpoint = exp(1i*angle_all); % all y' which on the endpoint
        re_yendpoint = real(y_endpoint);
        im_yendpoint = imag(y_endpoint);
        y_end = [re_yendpoint;im_yendpoint];

        V_V(1,:) = V.'*y_bold(:,1);
        index_plus = find(V_V(1,:)>0|V_V(1,:)==0);
        index_minus = find(V_V(1,:)<0);

        w_m(:,1) = sum(V(:,index_plus),2)-sum(V(:,index_minus),2); % sum by clown;

        for h = 2:(2*N(t)+2)
            V_V(h,:) = V.'*y_bold(:,h);
            index_n2p = find(sign(V_V(h,:))>0 & sign(V_V(h-1,:))<0);  %negative to positive
            index_p2n = find(sign(V_V(h,:))<0 & sign(V_V(h-1,:))>0);  %positive to negative

            w_m(:,h) = w_m(:,h-1)+ sum(V(:,index_n2p),2)-sum(V(:,index_p2n),2);
        end
        for h = 1:(2*N(t)+2)
            g_y1 = w_m(:,h).'*y_bold(:,h);
            g_y2 = w_m(:,h).'*y_end(:,h);
            if g_y1>=g_y2
                y_star(:,h) = y_bold(:,h);
            else
                y_star(:,h) = y_end(:,h);
            end
            g_ystar(h) = w_m(:,h).'*y_star(:,h);
        end
        index_y = find(max(g_ystar));
        y_final = y_star(:,index_y);
        for n = 1:N(t)
            x_star(n) = sign(V(:,1).'*y_final*V(:,n+1).'*y_final); % V is a matrix with dimension:2 x (N+1)
        end
        theta = acos(x_star);
        w_APX = exp(1i*theta);
        w_APX = w_APX.';
        y_square_APX(t) = [1;w_APX]'*R1*[1;w_APX];



        %%%Manopt
        manifold = complexcirclefactory(N(t)+1);
        problem.M = manifold;
        problem.cost = @(w) -w'*R*w; %
        problem.grad = @(w) manifold.egrad2rgrad(w,-2*R*w);
        %problem.grad = @(w) manifold.egrad2rgrad(w,-R.'*conj(w));%The mapping used here is a first-order orthogonal projection. detailed in "https://www.manopt.org/tutorial.html"
        %[w,wcost,info,options] = steepestdescent(problem); % at a random point on the manifold
        %[w,wcost,info,options] = conjugategradient(problem);
        [w,wcost,info,options] = barzilaiborwein(problem);  %Barzilai Borwein
        w = w./w(N(t)+1);
        w_2 = discretization1Bit(w);
        y_square_Man(t) = w_2'*R*w_2;
        w_Man = w_2(1:N(t));



        %%%SDR
                f_tmp = 0;
                r_tmp = zeros(N(t),1);   % r=[];The initialization of the three "tmp" variables cannot be placed outside the loop. They should be initialized whenever the value of "N" changes
                w_tmp = zeros(N(t),1);
        
                %CVX
                for k=1:count
                    r = (randn(N(t)+1,1)+1i*randn(N(t)+1,1)).*sqrt(1/2);   % (N,1)
                    cvx_begin
                    variable V(N(t)+1,N(t)+1) symmetric semidefinite   %
                    maximize( real(trace(R*V)))
                    subject to
                    diag(V) == 1;
                    cvx_end
        
                    [U,Sigma] = eig(V);
                    w = U*Sigma^(1/2)*r;   % (N*1)
                    f = w'*R*w;       %The number of random trials is given by 'count'. Obtain the beamforming vector w and Gaussian random vector r that correspond to the maximum value of f.
                    if f>f_tmp
                        f_tmp = max(f,f_tmp);
                        r_tmp = r;
                        w_tmp = w;     
                    end
                end
                w_tmp = w_tmp./w_tmp(N(t)+1);
                theta_opt = angle(w_tmp);
                w_sdr = exp(1i*theta_opt);
                w_SDR = discretization1Bit(w_sdr);
               % W = diag(w_SDR);
                y_square_SDR(t) = w_SDR'*R*w_SDR;
        
        %%%DaS
        [Z,eigs_R] = eig(R,'vector');
        [eigs_R,index] = sort(eigs_R,'descend');
        Z = Z(:,index);  %Arrange the vectors in Z according to their corresponding eigenvalues.
        Z_S = Z(:,1);      %The eigenvector corresponding to the maximum eigenvalue
        %%
        Phi_i = angle(Z_S);
        Phi_i(Phi_i>-pi & Phi_i<-pi/2) = Phi_i(Phi_i>-pi & Phi_i<-pi/2)+2*pi;
        I_1 = find(Phi_i>-pi/2 & Phi_i<pi/2);
        I_2 = find(Phi_i>pi/2 & Phi_i<3*pi/2);
        Phi_i_hat = Phi_i;
        Phi_i_hat(I_2) = Phi_i_hat(I_2)-pi;
        [Phi_i_hat,i_index] = sort(Phi_i_hat); %In ascending order
        %%  with Los
        W_hat = ones(N(t)+1,N(t)+1);
        for k = 1:N(t)+1
            W_hat(k,i_index(k+1:N(t)+1)) = -1;
        end
        for j = 1:N(t)+1
            W_hat(j,I_2) = -1*W_hat(j,I_2);
        end
        W = W_hat;
        W = W.';%Each w in W should be a column vector
        for e = 1:N(t)+1
            W(:,e) = W(:,e)./W(N(t)+1,e);%normalization
        end
        Y = W'*R*W;
        Y = diag(Y);%Take the elements on the diagonal
        [y_square_DaS(t),w_index] = max(Y);
        w_1 = W(:,w_index);
        w_DaS = w_1(1:N(t));
    end
    y_APX(m,:) = sqrt(y_square_APX);
    y_Man(m,:) = sqrt(y_square_Man);
    y_SDR(m,:) = sqrt(y_square_SDR);
    y_DaS(m,:) = sqrt(y_square_DaS);
    X = sprintf('The loop have completed %d times.',m);
    disp(X);
end
yAPX = mean(y_APX,1); 
yMan = mean(y_Man,1);
ySDR = mean(y_SDR,1); 
yDaS = mean(y_DaS,1); 

SNR_APX = 10*log10(1+yAPX.^2./sigma2);
SNR_Man = 10*log10(1+yMan.^2./sigma2); 
SNR_SDR = 10*log10(1+ySDR.^2./sigma2);
SNR_DaS = 10*log10(1+yDaS.^2./sigma2);

figure
 hold on
 plot(N,SNR_APX,'-kd','LineWidth',2);
 hold on 
 plot(N,SNR_Man,'-.bo','LineWidth',2);
 hold on
 plot(N,SNR_SDR,'--g*','LineWidth',2);
 hold on 
 plot(N,SNR_DaS,'-rs','LineWidth',2);
 hold on
 xlabel('N','Interpreter','latex','Fontsize',16);
 ylabel('SNR(dB)','Interpreter','latex','Fontsize',16);
 legend('Solved by APX','Solved by Manopt','Solved by SDR','Solved by Proposed DaS','Interpreter','latex','Fontsize',12);
 hold on ;
 grid on;
 box on;
