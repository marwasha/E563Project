% Simulate 3 bus Network of Synchronous Machines
opts = odeset('RelTol',1e-6,'AbsTol',1e-8, 'Vectorized', 'off');
traj = [];
times = [];
x0_NL = [0.3358   -0.0000    0.3029    0.0000   -0.0460   -0.0000];
p = [0.5,1,1.2,0.5,0.5,0.5];
u = [0.4,0.4,-0.8,1,1,1];
uF = [0.4,0.0,-0.8,1,1,1];

simtime = 1000;
faultTime = .04;
[time, X_NL] = ode23t(@(t, x) grid_model(t, x,uF,p), [0 faultTime], x0_NL,opts);
x0_NL = X_NL(end,:);
times = [times;time];
traj = [traj;X_NL];
[time, X_NL] = ode23t(@(t, x) grid_model(t, x,u,p), [faultTime simtime], x0_NL,opts);
times = [times;time];
traj = [traj;X_NL];
figure

plot(traj(:,1),traj(:,3))

dyn = sym_gen_2;
A = full(dyn.A(x0_NL,u,p));
B = full(dyn.B(x0_NL,u,p));