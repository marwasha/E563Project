function out = reachabilityAnaylsis(X0, U0, P0, param)

dynParam = nonlinearSys('dyn', @grid_model_CORA, 12, 6);

paramsFault.tFinal = param.fault.tFinal;
paramsFault.R0 = zonotope([[X0.init; P0.init], diag([X0.dist, P0.dist])]);
paramsFault.U = zonotope([U0.fault.init, diag(U0.fault.dist)]);

out.fault.R0 = paramsFault.R0; 

optionsFault.timeStep= paramsFault.tFinal/param.fault.steps;
optionsFault.taylorTerms=4;
optionsFault.intermediateOrder = 4;
optionsFault.zonotopeOrder=10;
optionsFault.tensorOrder = 2;
optionsFault.alg = 'lin';

tic
RcontFault = reach(dynParam, paramsFault, optionsFault);
tComp = toc;
disp(['computation time of reachable set during fault: ',num2str(tComp)]);

simOpt.points = param.points;        % number of initial points
simOpt.fracVert = 0.5;     % fraction of vertices initial set
simOpt.fracInpVert = 0.5;  % fraction of vertices input set
simOpt.inpChanges = 0;     % changes of input over time horizon

simResFault = simulateRandom(dynParam,paramsFault,simOpt);

out.fault.reach = RcontFault;
out.fault.sim = simResFault;

paramsPost.tFinal = param.post.tFinal;
paramsPost.R0 = RcontFault.timePoint.set{end}; %Initialize at final set of last
paramsPost.U = zonotope([U0.post.init, diag(U0.post.dist)]);

out.post.R0 = paramsPost.R0;

optionsPost.timeStep= paramsPost.tFinal/param.post.steps;
optionsPost.taylorTerms=4;
optionsPost.intermediateOrder = 4;
optionsPost.zonotopeOrder=12;
optionsPost.tensorOrder = 2;
optionsPost.alg = 'lin';

tic
RcontPost = reach(dynParam, paramsPost, optionsPost);
tComp = toc;
disp(['computation time of reachable set after fault: ',num2str(tComp)]);

simResPost = simulateRandom(dynParam,paramsPost,simOpt);

out.post.reach = RcontPost;
out.post.sim = simResPost;

end