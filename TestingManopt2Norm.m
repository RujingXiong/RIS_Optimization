clear all;close all;clc;

%% Generate random data.
numBits = 1;

stepAngle = 2*pi/2^numBits;

n = 100;

m = 10;

p = 2;

numTrials = 20000;


%%
manifold = complexcirclefactory(n);
problem.M = manifold;

costContinuous = zeros(numTrials, 1);
costRounded = zeros(numTrials, 1);
costLifting = zeros(numTrials, 1);

f = waitbar(0,'please wait...');

for iTrial = 1:1:numTrials

    waitbar(iTrial/numTrials,f,'please wait...');

    A = randn(m,n) + 1i*randn(m,n);

    A_A = A'*A;

    % Define the problem cost function and its Euclidean gradient.
    problem.cost  = @(x) -x'*(A_A*x);
    problem.egrad = @(x) -2*A_A*x;      % notice the 'e' in 'egrad' for Euclidean

    [x, xcost, info, options] = trustregions(problem);
    
    xAlternatingRounded = exp(1i*stepAngle.*round(angle(x)./stepAngle));
    
    [xLifting, costList2] = Lift_max_Ax_p(A, p, xAlternatingRounded, numBits);
    

    costContinuous(iTrial) = norm(A*x, p);
    costRounded(iTrial) = norm(A*xAlternatingRounded, p);
    costLifting(iTrial) = norm(A*xLifting, p);

end

close(f)

mean(costLifting - costRounded)

%%
nameWorkspace = strcat('WS_2Norm-', datestr(now,'mmmm-dd-HH-MM-SS'), '.mat');

save(nameWorkspace)

%%
relativeLiftingGain = abs(costLifting-costRounded)./abs(costContinuous-costRounded);

%%
figure(1)
histogram(costLifting-costRounded, 'Normalization','count')
grid on
xlim([-0.3, 9])

xlabel('$|| A e^{j \mathbf{\Omega}_{lift} } ||_2 - || A e^{j \mathbf{\Omega}_{round} } ||_2$','Interpreter','latex','FontSize',12)
ylabel('Count','FontSize',12)

% exportgraphics(gcf, 'Lifting_2Norm_20000.pdf')

%%
figure(2)
histogram(relativeLiftingGain, 'Normalization','count')
grid on
xlim([0, 0.50])
xline(median(relativeLiftingGain),'-r','LineWidth',2)


xlabel('$\frac{|| A e^{j \mathbf{\Omega}_{lifted} } ||_2 - || A e^{j \mathbf{\Omega}_{rounded} } ||_2}{|| A e^{j \mathbf{\Omega}_{unrounded}  } ||_2 - || A e^{j \mathbf{\Omega}_{rounded} } ||_2}$','Interpreter','latex')
ylabel('Count')

%matlab2tikz('file1.tex');

% exportgraphics(gcf, 'Lifting_2Norm_20000_Percent.pdf')