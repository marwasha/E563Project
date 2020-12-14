addpath(genpath('ToolBoxes\CORA_2020'))
clear all

%% Set up
%[0 0 0.4691 0 0.1236 0]
X0.init = [0.0000; 0.0000; 0.4691; 0.0000; 0.1236; 0.0000];
X0.dist = [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000];
P0.init = [0.5000; 1.0000; 1.2000; 0.5000; 0.5000; 0.5000];
P0.dist = [0.0500, 0.1000, 0.1200, 0.0500, 0.0500, 0.0500]/4;

U0.fault.init = [-0.5; 0.0; -0.3; 1.0; 1.0; 1.0];
U0.fault.dist = [ 0.0000, 0.000, 0.0000, 0.0000, 0.0000, 0.0000];
param.fault.tFinal = 0.4;
param.fault.steps = 100;

U0.post.init = [-0.5; 0.8; -0.3; 1.0; 1.0; 1.0];
U0.post.dist = [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000];
param.post.tFinal = .0000000001;
param.post.steps = 10;

param.points = 1;

%% Run

res = reachabilityAnaylsis(X0, U0, P0, param);

%% Plot

dims = {[1,2],[2,3],[1,3]};

for i = 1:length(dims)
    
    figure; hold on; box on;
    projDims = dims{i};

    % plot reachable sets
    hanFault = plot(res.fault.reach,projDims,'FaceColor',[1 0 0],'EdgeColor','none');
    hanPost = plot(res.post.reach,projDims,'FaceColor',[0 0 1],'EdgeColor','none', 'FaceAlpha', 1);
    
    % plot initial set
    plot(res.fault.R0,projDims,'w','Filled',true,'EdgeColor','k');
  
    % plot simulation results
    plot(res.fault.sim,projDims);
    plot(res.post.sim,projDims);
    %label plot
    xlabel(['x_{',num2str(projDims(1)),'}']);

    ylabel(['x_{',num2str(projDims(2)),'}']);
    legend([hanFault, hanPost],'fault', 'post');
end