clc; clear all;

%% Some parameters
params.gridSize = 10;
params.a = 0.1;             % Coefficient a of the objective function
params.eps = 0.1;           % Small dampening added to the laplacian before
                            % computing the eigenvalues
params.minWeight = 0.2;     % Constraint for weights

%% Initialization

w0 = ones(params.gridSize*(params.gridSize-1)*2,1); % Initial guess, same 
                                                    % as the initial state
                                                    % of the graph

% objFuncHandle = @objectiveFunctionWithGradient;
objFuncHandle = @altObjectiveFunctionWithGradient;

disp(['Initial Objective (old): ' ...
    num2str(objectiveFunctionWithGradient(w0,params))])

disp(['Initial Objective (alt): ' ...
    num2str(altObjectiveFunctionWithGradient(w0,params))])

% No linear inequality constraints
A = [];
b = [];

% Weight sum constraint implemented as a linear constraint
Aeq = ones(1,length(w0));
beq = sum(w0);

% It is possible to provide lower and upper bounds for the weights
lb = params.minWeight*ones(length(w0),1);
ub = inf*ones(length(w0),1);

% No nonlinear constraints
nonlcon = [];

% Use options for utilizing the known gradient
options = optimoptions('fmincon','SpecifyObjectiveGradient',true);

w = fmincon(@(w)objFuncHandle(w,params),w0,A,b,Aeq,beq,...
    lb,ub,nonlcon,options);

disp(['Final Objective (old): ' ...
    num2str(objectiveFunctionWithGradient(w,params))])

disp(['Final Objective (alt): ' ...
    num2str(altObjectiveFunctionWithGradient(w,params))])

%% Post Processing for Generating Histograms

[A_init,D_init,L_init] = generateMatricesFromWeights(w0,params.gridSize);
[~,diag_lambda_init] = eig(L_init + params.eps*eye(length(L_init)));
lambda_init = diag(diag_lambda_init);

[A_final,D_final,L_final] = generateMatricesFromWeights(w,params.gridSize);
[~,diag_lambda_final] = eig(L_final + params.eps*eye(length(L_final)));
lambda_final = diag(diag_lambda_final);

figure()
subplot(2,1,1)
hist_init = histogram(lambda_init,10);
maxFreq = max(hist_init.BinCounts);
ylim([0,maxFreq])
title("Initial spectrum")
xlabel("Eigenvalue")
ylabel("Frequency")
subplot(2,1,2)
hist_final = histogram(lambda_final,10);
ylim([0,maxFreq])
title("Final spectrum")
xlabel("Eigenvalue")
ylabel("Frequency")
