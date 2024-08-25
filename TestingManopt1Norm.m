clear all;close all;clc;

%% Generate random data.
numBits = 1;

stepAngle = 2*pi/2^numBits;

n = 100;

m = 10;

p = 1;

numTrials = 20000;

%% Create the problem structure.

manifold = complexcirclefactory(n);
problem.M = manifold;

costContinuous = zeros(numTrials, 1);
costRounded = zeros(numTrials, 1);
costLifting = zeros(numTrials, 1);

f = waitbar(0,'please wait...');

for iTrial = 1:1:numTrials

    waitbar(iTrial/numTrials,f,'please wait...');

    A = randn(m,n) + 1i*randn(m,n);
    
    % Define the problem cost function and its Euclidean gradient.
    problem.cost  = @(x) -norm(A*x, 1);
    problem.egrad = @(x) -A'*sign(A*x);      % notice the 'e' in 'egrad' for Euclidean
 
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
nameWorkspace = strcat('WS_1Norm-', datestr(now,'mmmm-dd-HH-MM-SS'), '.mat');

save(nameWorkspace)


%%
relativeLiftingGain = abs(costLifting-costRounded)./abs(costContinuous-costRounded);


%%
figure(1)
histogram(costLifting-costRounded, 'Normalization','count')
xlim([-0.5, 30])

xlabel('$|| A e^{j \mathbf{\Omega}_{lift} } ||_1 - || A e^{j \mathbf{\Omega}_{round} } ||_1$','Interpreter','latex','FontSize',12)
ylabel('Count','FontSize',12)
grid on

% exportgraphics(gcf, 'Lifting_1Norm_20000.pdf')

%%
figure(2)
histogram(relativeLiftingGain, 'Normalization','count')
grid on
xlim([0, 0.50])
xline(median(relativeLiftingGain),'-r','LineWidth',2)

xlabel('$\frac{|| A e^{j \mathbf{\Omega}_{lifted} } ||_1 - || A e^{j \mathbf{\Omega}_{rounded} } ||_1}{|| A e^{j \mathbf{\Omega}_{unrounded} } ||_1 - || A e^{j \mathbf{\Omega}_{rounded} } ||_1}$','Interpreter','latex')
ylabel('Count')

% exportgraphics(gcf, 'Lifting_1Norm_20000_Percent.pdf')

matlab2tikz('file1.tex');