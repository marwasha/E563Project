% Simulate 3 bus Network of Synchronous Machines
opts = odeset('RelTol',1e-6,'AbsTol',1e-8, 'Vectorized', 'off');
dyn = Lyap_sym_gen();
traj = [];
times = [];
lyapunov = [];
lyapunovdx = [];
SEP = [0.3358 0 0.3029 0 -0.0460 0];
x0_NL = [0.3358 0 0.3029 0 -0.0460 0];
p = [0.5,1,1.2,.1,.1,.1];
u = [0.4,0.4,-0.8,1,1,1];
uF = [0.4,0.0,-0.8,1,1,1];

simtime = 101;
faultTime = .00000001;
[time, X_NL] = ode23t(@(t, x) grid_model(t, x,uF,p), [0 faultTime], x0_NL,opts);
x0_NL = X_NL(end,:);
times = [times;time];
stepsNL = length(time);
traj = [traj;X_NL];
[time, X_NL] = ode23t(@(t, x) grid_model(t, x,u,p), [faultTime simtime], x0_NL,opts);
times = [times;time];
traj = [traj;X_NL];

for i = 1:length(traj)
    if i < stepsNL
        uL = uF;
    else
        uL = u;
    end
    [Lyapunov,PE,KE] = lyapunov_2(traj(i,:),SEP,uL(1:3));
    Lyapunovdx = full(dyn.dvdx(traj(i,:),SEP,uL(1:3)))*grid_model(1,traj(i,:),uL,p);
    lyapunov = [lyapunov, Lyapunov];
    lyapunovdx = [lyapunovdx, Lyapunovdx];
end

figure
plot(traj(:,1)-traj(:,3), traj(:,1)-traj(:,5))
title("traj")
figure
plot(lyapunov)
title("lyap")
figure
plot(lyapunovdx)
title("lyapdv")
dyn = sym_gen_2;
A = full(dyn.A(SEP,u,p));
B = full(dyn.B(SEP,u,p));

%% 

P = sdpvar(6,6);
gamma = sdpvar(1,1);
Obj = [];
eps = 10^-15;

Consts = [P>= eps*eye(6), A'*P + P*A <= -eps*eye(6)];

ops = sdpsettings('solver', 'sedumi');

diagnostics = solvesdp(Consts, Obj, ops);

if diagnostics.problem == 0
 disp('Feasible: The augmented system is stable')
 Pval = double(P)
 eig(Pval)
else
 disp('No quadratic common Lyapunov function') 
 disp('May or may not be stable')
end


