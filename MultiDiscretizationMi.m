clear all;clc;%close all;
sigma2 = 1;  %noise
T = 1;  %number of Tx
N = 200:200:1000; %number of RIS elements
loop = 100;
%%

numBits = 5;

%%
angleWeightRange = 2*pi/2^numBits;
y_my = zeros(loop,length(N)); % my method
y_man = zeros(loop,length(N)); % manopt
for h = 1:loop
    for t = 1:length(N)
        H = (randn(N(t),T)+1i*randn(N(t),1)).*0.5;     %AP-RIS
        Hr = (randn(N(t),1)+1i*randn(N(t),1));%RIS-USER
        Phi = diag(Hr)*H;
        R = Phi*Phi';
        [Z,eigs_R] = eig(R,'vector');
        [eigs_R,index_z] = sort(eigs_R,'descend');
        Z = Z(:,index_z);  %The corresponding eigenvectors are sorted in descending order too.
        z = Z(:,1);      %Taking the eigenvector corresponding to the largest eigenvalue
        %%
        %z = randn(N,1)+1i.*randn(N,1);
        %z = z./norm(z);

        %%
        thetaZ = mod(angle(z),2*pi);

        thetaZNormalized = mod(thetaZ,angleWeightRange);

        [~, indexThetaZ] = sort(thetaZNormalized,'ascend');

        %%
        thetaZShift = round((thetaZ-thetaZNormalized)/angleWeightRange);
        thetaZShifting = repmat(thetaZShift.',2^numBits*N(t),1);

        %%
        angleWeightNormalized = zeros(N(t)*2^numBits,N(t));

        for iPart = 1:1:2^numBits
            angleWeightNormalized((iPart-1)*N(t)+indexThetaZ, indexThetaZ) = (iPart-1)*ones(N(t),N(t)) + triu(ones(N(t),N(t))).';
        end

        angleWeightNormalized = angleWeightNormalized - thetaZShifting;

        angleWeightNormalized = mod(angleWeightNormalized,2^numBits);

        %%
        weightReduced = exp(1i*angleWeightRange.*angleWeightNormalized);

        [valueOpt, indexOpt] = max(abs(weightReduced*z));

        %%
        weightOpt = weightReduced(indexOpt,:);
        weightOpt1 = weightOpt';
        y_my(h,t) = weightOpt1'*R*weightOpt1;
%       y_my(h,t) = abs(weightOpt*z);
        %% Manopt
        manifold = complexcirclefactory(N(t));
        problem.M = manifold;
        problem.cost = @(w) -w'*R*w; %
        problem.grad = @(w) manifold.egrad2rgrad(w,-2*R*w);
        %problem.grad = @(w) manifold.egrad2rgrad(w,-R.'*conj(w));%
        [w,wcost,info,options] = steepestdescent(problem); % at a random point on the manifold
        %[w,wcost,info,options] = conjugategradient(problem);%
        %[w,wcost,info,options] = barzilaiborwein(problem);  %Barzilai Borwein
        % w = Discretization(w);
        %R = z*z';
%         y_man(h,t) =abs(w'*z);
        y_man(h,t) =w'*R*w;
    end
    X = sprintf('The loop have completed %d times.',h);
    disp(X);
end
y_meanMy = mean(y_my,1);
SNR_my = 10*log10(1+y_meanMy./sigma2);

y_meanMan = mean(y_man,1);
SNR_man = 10*log10(1+y_meanMan./sigma2);

figure
hold on
plot(N,SNR_man,'-rs','LineWidth',2);
hold on
plot(N,SNR_my,'-bs','LineWidth',2);
xlabel('N(t)','Interpreter','latex','Fontsize',15);
ylabel('SNR(dB)','Interpreter','latex','Fontsize',15);
legend('Solved by Manopt','Solved by proposed method');
hold on ;
grid on;
hold on;
box on;
title('Gains on Signal power');
toc

%%
figure
hold on
plot(N,SNR1,'-bs','LineWidth',1);
hold on
plot(N,SNR2,'-cd','LineWidth',1);
hold on
plot(N,SNR3,'-r*','LineWidth',1);
hold on
plot(N,SNR4,'-cx','LineWidth',1);
hold on
plot(N,SNR5,'--ks','LineWidth',1);
hold on
plot(N,SNR6,'--bs','LineWidth',1);
hold on
plot(N,SNR7,'--rd','LineWidth',1);
hold on
plot(N,SNR8,'--co','LineWidth',1);
xlabel('N','Interpreter','latex','Fontsize',18);
ylabel('SNR(dB)','Interpreter','latex','Fontsize',18);
legend('1-bit','2-bit','3-bit','4-bit','5-bit','6-bit','7-bit','8-bit','Interpreter','latex','Fontsize',16);
box on;
hold on
grid on;