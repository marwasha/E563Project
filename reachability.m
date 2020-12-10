addpath(genpath('ToolBoxes\CORA_2020'))
clear all
% Parameter ---------------------------------------------------------------

params.tFinal = 10;
params.R0 = zonotope([[0; 0; 0.4691; 0; 0.1236; 0],0.01*eye(6)]);
params.U = zonotope([[-0.5;0.8;-0.3;1;1;1],0.01*eye(6)]);

% Reachability Settings ---------------------------------------------------

options.timeStep=0.5;
options.taylorTerms=4;
options.intermediateOrder = 4;
options.zonotopeOrder=10;
options.tensorOrder = 2;
options.alg = 'lin';

% System Dynamics ---------------------------------------------------------

% tank system with uncertain parameters
optionsParam = options;
optionsParam.paramInt = interval([0.5;1;1.2;0;0;0]-.01,[0.5;1;1.2;0;0;0]+.01);

dynParam = nonlinParamSys('dyn', @grid_model_CORA, 6, 6, 6, 'constParam');

% Reachability Analysis ---------------------------------------------------    

% compute reachable set of tank system with uncertain parameters
tic
RcontParam = reach(dynParam, params, optionsParam);
tComp = toc;
disp(['computation time of reachable set with uncertain parameters: ',num2str(tComp)]);