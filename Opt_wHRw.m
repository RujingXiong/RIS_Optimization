function[w_opt, valueOpt] = Opt_wHRw(R,B)
%% Author: Rujing Xiong, Tiebin Mi
 % CreatTime: 2022.09.26
 % Complete:2022.09.26
 % Modified: 2022.09.30 modified for multi-bit
 % E-mail: Rujing@hust.edu.cn
 % Description: DiscreteOptimization for 2-bit RIS, multi-bit
 % Based on traditional channel model
 % R is rank-one
%% check input argument
if ~exist('R','var')
    error('The Z is empty.')
end

if ~exist('B','var')
    warning('The numBits is empty. Set numBits=4 by default.')
    B = 4;
end
%%
[Z,eigs_R] = eig(R,'vector');
[~,index_z] = sort(eigs_R,'descend');
Z = Z(:,index_z);  %The corresponding eigenvectors are sorted in descending order too.
z = Z(:,1);      %Taking the igenvector corresponding to the largest eigenvalue

N = length(z);
U = zeros(N,2^B*N);

theta_ori = angle(z);  % the angle of eigenvector
theta_temp = angle(z);
theta_ori = mod(theta_ori,2*pi);% into(0,2*pi]


theta_tilde = repmat(theta_ori,1,2^B);
partition_maxrix = 0:2*pi/2^B:(2*pi-2*pi/2^B);

theta_tilde = theta_tilde + partition_maxrix;%compute theta_tilde

theta_tilde = mod(theta_tilde,2*pi);% devote theta_tilde (0,2*pi]
theta_tilde = sort(theta_tilde,2); % sort each row, select the smallest candidate from the four

theta = theta_tilde(:,1);% take the first column as initialization

[theta,index_t] = sort(theta); % sort theta_tilde. Note! need to record the first N(t) elements' index. Because it include the solution index imformation

U_ = repmat(theta,1,2^B*N); % repeat the row of theta_tilde once, the column 4*N(t) times to obtain U.
%theta_hat = sort(reshape(theta_tilde,2^B*N(t),1))+pi/4;


triu_temp = triu( (2*pi/2^B)*ones(N) );% a upper triangular matrix which element is 2*pi/2^B
triu_temp = repmat(triu_temp,1,2^B) ;% repeat the upper triangular 2^B times by row
U_partition_temp = kron(( 0:1:(2^B-1) ).*(2*pi/2^B),ones(N) ); %calculate the kronecker product;
U_partition = U_partition_temp + triu_temp;
U_hat = U_partition + U_;

for i = 1:N
    U(index_t(i),:) = U_hat(i,:);%Thereby, set U is completed.
end

temp = repmat(theta_temp,1,2^B*N);
Tau = U-temp;% Tau is the matrix that each colum of U minus theta_ori; each column of Tau is the angles of a candidate w.

Tau = cos(Tau)+1i*sin(Tau);%above tau in Tau is phase.
%Tau = conj(Tau);
Y = Tau'*R*Tau;
Y = diag(Y);% get the diagonal elements
[~,Tau_index] = max(Y);
w_opt = Tau(:,Tau_index);
w_opt = conj(w_opt);
valueOpt =w_opt'*R*w_opt;