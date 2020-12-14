function simulator(p)

theta_eq = asin(p(1)/p(2));
SEP = theta_eq;
UEP = pi-theta_eq;
UEP_w = 0;
tspan = 0:1:15;
traj1 = [];
epsilon = 1e-4;
tspan = -15:1:0;
func = @(t,x) [-x(2);-1*(P_m - P*sin(x(1))-D*x(2))];
Initial = kron([UEP,0],ones(10,1));
Initial_2 = kron([UEP,0],ones(10,1));
Initial = Initial + [-1,1].*abs(epsilon*randn(10,2));
Initial_2 = Initial_2 + [1,-1].*abs(epsilon*randn(10,2));
Initial = [Initial;Initial_2];

for k = 1:1:length(Initial)
    X_0 = Initial(k,:);
    traj_2 = [];
    for i = 1:1:length(tspan)-1
        ts = [tspan(i) tspan(i+1)];
        [t,x] = ode45(@(t,x)single_machine_model_backwards(t,x,p), ts, X_0);
        X_0 = x(end,:);
        traj_2 = [traj_2;x];
    end
    plot(traj_2(:,1),traj_2(:,2),'b')
    hold on
end

grid on
axis equal
xlim([-2 5])
ylim([-4 3])
xlabel('\theta (rad)')
ylabel('\omega (rad/s)')
hold on
plot(SEP,0,'rx')
legend('Region of Attraction','Location','northeastoutside')
hold on

X_0 = [1;-1];
[t,x] = ode45(@(t,x)single_machine_model(t,x,p), [0 50], X_0);
plot(x(:,1),x(:,2))

% Check for P satisfying Lya
dyn = sym_gen_single();
A = full(dyn.A([0;-2],p));

P = sdpvar(2,2);
Obj = [];
eps = 10^-8;

Consts = [P>= eps*eye(2), A'*P + P*A <= -eps*eye(2)];

ops = sdpsettings('solver', 'sedumi');

diagnostics = solvesdp(Consts, Obj, ops);

if diagnostics.problem == 0
 disp('Feasible: The augmented system is stable')
 Pval = double(P)
else
 disp('No quadratic common Lyapunov function') 
 disp('May or may not be stable')
end

eig(Pval)
end

