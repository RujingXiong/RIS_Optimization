clear all;clc;

%%
numBits = 2;

stepAngle = 2*pi/2^numBits;

n = 100;

m = 30;

A = randn(m,n) + 1i*randn(m,n);


load('A.mat');

%%
w_k = exp(1i*2*pi*rand(1))*ones(n,1);

w_k1 = exp(1i*2*pi*rand(n,1));


iIter = 0;

costValue = [];

while abs(norm(A*w_k1,1) - norm(A*w_k,1))/norm(A*w_k1,1) > 1e-6

    iIter = iIter + 1;

    w_k = w_k1;

    z_k = exp(1i*angle(A*w_k));
    w_k1 = exp(1i*angle(A'*z_k));


%     z_k = A*w_k./norm(A*w_k);
%     w_k1 = exp(1i*angle(A'*z_k));

    costValue(iIter) = norm(A*w_k1,1);

end

%%
w_k1_Quant = exp(1i*stepAngle.*round(angle(w_k1)./stepAngle));

iIter = iIter+1;

costValue(iIter) = norm(A*w_k1_Quant,1);

%%
w_k = exp(1i*2*pi*rand(1))*ones(n,1);

w_k1 = w_k1_Quant;
k = 1;
while k<15
%while 1e10*abs(norm(A*w_k1,1) - norm(A*w_k,1))/norm(A*w_k1,1) > 1e-10

    iIter = iIter + 1;

    w_k = w_k1;
    
    z_k = exp(1i*angle(A*w_k));
    w_k1 = Opt_absWZ(A'*z_k, numBits);

%     z_k = exp(1i*2*pi*rand(1))*exp(1i*angle(A*w_k));
%     w_k1 = Opt_absWZ(A'*z_k, numBits);

    norm(A*w_k1,1)

    costValue(iIter) = norm(A*w_k1,1);
k = k+1;
end


%%
figure
plot(costValue, 'Marker', 'square', 'LineWidth', 1, 'LineStyle', '-', 'Color', [0 0 1]);
grid on
xlim([1, length(costValue)])
% ylim([400, 480])

xlabel('$k$','Interpreter','latex')
ylabel('$|| A e^{j \mathbf{\Omega} } ||_1$','Interpreter','latex')


%%
% matlab2tikz('file2.tex');