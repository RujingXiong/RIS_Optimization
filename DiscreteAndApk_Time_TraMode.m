%% Author: Rujing Xiong
 % Time: 2022.07.16
 % Complete:2022.07.29
 % Modified:2022.07.30
 % E-mail: Rujing@hust.edu.cn
 % Description: Reappear the DiscreteOptimization in paper 'Configuring intelligent reflecting surface with performance guarantees: Optimal beamforming'
 % And execution-time Comparison
 %%  SNR with N
 clc
 clear 
 tic
 %sigma2 = 1;  %noise
 T = 1;  %number of Tx
 N = 100;
 loop = 100;
 for h = 1:loop
     for t = 1:length(N)
         w_m = zeros(2,2*N(t)+2);
         y_star= zeros(2,2*N(t)+2);
         x_star = zeros(1,N(t));
         g_ystar = zeros(1,2*N(t)+2);
         V_V = zeros(1,N(t)+1);
         H = (randn(N(t),T)+1i*randn(N(t),1)).*2;     %AP-RIS
         Hr = (randn(N(t),1)+1i*randn(N(t),1)).*1;%RIS-USER
         % Hr = [3*ones(N(t)/3,1); 1*ones(N(t)/3,1);5*ones(N(t)/3,1)];%RIS-USER

         G = (randn(1,T)+1i*randn(1,T)).*sqrt(1/2).*1;  %LOS

         %         Phi = diag(Hr)*H;
         h_m = [diag(Hr)*H;G]';
         Re_hm = real(h_m);
         Im_hm = imag(h_m);
         V = [Re_hm;Im_hm];
         v_0 = V(:,1);

         R = h_m'*h_m;

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

         for m = 2:(2*N(t)+2)
             V_V(m,:) = V.'*y_bold(:,m);  
             index_n2p = find(sign(V_V(m,:))>0 & sign(V_V(m-1,:))<0);  %negative to positive
             index_p2n = find(sign(V_V(m,:))<0 & sign(V_V(m-1,:))>0);  %positive to negative

             w_m(:,m) = w_m(:,m-1)+ sum(V(:,index_n2p),2)-sum(V(:,index_p2n),2);
         end
         for m = 1:(2*N(t)+2)
             g_y1 = w_m(:,m).'*y_bold(:,m);
             g_y2 = w_m(:,m).'*y_end(:,m);
             if g_y1>=g_y2
                 y_star(:,m) = y_bold(:,m);
             else
                 y_star(:,m) = y_end(:,m);
             end
             g_ystar(m) = w_m(:,m).'*y_star(:,m);
         end
         index_y = find(max(g_ystar));
         y_final = y_star(:,index_y);

         for n = 1:N(t)
             x_star(n) = sign(V(:,1).'*y_final*V(:,n+1).'*y_final); % V is a matrix with dimension:2 x (N+1)
         end
         theta = acos(x_star);
         w = exp(1i*theta);
         w = w.';
     end
     X = sprintf('The loop have completed %d times.',h);
     disp(X);
 end
toc
 %% Author: Rujing Xiong
 % Time: 2022.07.16
 % Complete:2022.07.29
 % Modified:2022.07.30
 % E-mail: Rujing@hust.edu.cn
 % Description: Reappear the DiscreteOptimization in paper 'Configuring intelligent reflecting surface with performance guarantees: Optimal beamforming'
 % And execution-time Comparison
 %%  SNR with N
 clc
 clear
 tic
 sigma2 = 1;  %noise
 T = 1;  %number of Tx
 N = 100;

 y_square = zeros(1,length(N));

 loop = 100;
 y_luo = zeros(loop,length(N)); % luo's method
 y_my = zeros(loop,length(N)); % my method

 for h = 1:loop
     for t = 1:length(N)
         w_m = zeros(2,2*N(t)+2);
         y_star= zeros(2,2*N(t)+2);
         x_star = zeros(1,N(t));
         g_ystar = zeros(1,2*N(t)+2);
         V_V = zeros(1,N(t)+1);
         H = (randn(N(t),T)+1i*randn(N(t),1)).*2;     %AP-RIS
         Hr = (randn(N(t),1)+1i*randn(N(t),1)).*1;%RIS-USER
         % Hr = [3*ones(N(t)/3,1); 1*ones(N(t)/3,1);5*ones(N(t)/3,1)];%RIS-USER

         G = (randn(1,T)+1i*randn(1,T)).*sqrt(1/2).*1;  %LOS

         %         Phi = diag(Hr)*H;
         h_m = [diag(Hr)*H;G]';
         Re_hm = real(h_m);
         Im_hm = imag(h_m);
         V = [Re_hm;Im_hm];
         v_0 = V(:,1);

         R = h_m'*h_m;
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
        W = W.';% every w in W is a column vector
        for e = 1:N(t)+1
            W(:,e) = W(:,e)./W(N(t)+1,e);%each w keeps the last element 1
        end
        Y = W'*R*W;
        Y = diag(Y);% get the diagonal elements
        [y_square(t),w_index] = max(Y);
        w_opt = W(:,w_index);
        w_my = w_opt(1:N(t));
     end
     X = sprintf('The loop have completed %d times.',h);
     disp(X);
 end
 toc