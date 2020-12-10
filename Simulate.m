% Simulate 3 bus Network of Synchronous Machines
opts = odeset('RelTol',1e-6,'AbsTol',1e-8, 'Vectorized', 'off');

%x0_NL = [0 0 0.4691 0 0.1236 0]; EQ when was wrt t1,t2,t3
x0_NL = [-.4691 -0.1236 0.3455 0 0 0]; 
uStab = [-0.5,0.8,-0.3,1,1,1];
uFault = [-0.5,0.0,-0.3,1,1,1];
p = [0.5,1,1.2,0,0,0];

t_step = .001;
t_fault = .004;
totalSimTime = 2;
steps = round(totalSimTime/t_step);
stepFault = round(t_fault/t_step);

traj = [];
times = [];

[time, X_NL] = ode23t(@(t, x) grid_model(t, x, uFault, p), [0 t_fault] , x0_NL,opts);
traj = [traj; X_NL];
times = [times; time];
[time, X_NL] = ode23t(@(t, x) grid_model(t, x, uStab, p), [t_fault totalSimTime],X_NL(end,:),opts);
traj = [traj; X_NL];
times = [times; time];

figure
plot(traj(:,1), traj(:,2))
figure
plot(traj(:,2), traj(:,3))
figure
plot(traj(:,1), traj(:,3))

dyn = sym_gen_2;
A = full(dyn.A(x0_NL, uStab,p));
B = full(dyn.B(x0_NL, uStab,p));