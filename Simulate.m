% Simulate 3 bus Network of Synchronous Machines
opts = odeset('RelTol',1e-6,'AbsTol',1e-8, 'Vectorized', 'off');
traj = [];
SEP = [0.3358 0 0.3029 0 -0.0460 0];
x0_NL = [0.3 0 0.30 0 0 0];
p = [0.5,1,1.2,0.1,0.1,0.1,0.4,0.4,-0.8,1,1,1];
lyapunov = [];
for i = 1:1:1000
    [time, X_NL] = ode23t(@(t, x) grid_model(t, x,p), 0.01:0.01:0.1, x0_NL,opts);
    x0_NL = X_NL(end,:);
    traj = [traj;X_NL];
    [Lyapunov,PE,KE] = lyapunov_2(x0_NL,SEP,p(7:9));
    lyapunov = [lyapunov,Lyapunov];
end
figure
plot(traj(:,1),traj(:,3))
figure
plot(lyapunov)
dyn = sym_gen_2;
A = full(dyn.A(x0_NL,p));