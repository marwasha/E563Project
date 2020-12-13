% Simulate 3 bus Network of Synchronous Machines
opts = odeset('RelTol',1e-6,'AbsTol',1e-8, 'Vectorized', 'off');
traj = [];
x0_NL = [0 0 0.4691 0 0.1236 0];
p = [0.5,1,1.2,0.5,0.5,0.5,0.4,0.4,-0.8,1,1,1];
for i = 1:1:1000
    [time, X_NL] = ode23t(@(t, x) grid_model(t, x,p), 0.01:0.01:0.4, x0_NL,opts);
    x0_NL = X_NL(end,:);
    traj = [traj;X_NL];
end
figure
plot(traj(:,1),traj(:,3))

dyn = sym_gen_2;
A = full(dyn.A(x0_NL,p));