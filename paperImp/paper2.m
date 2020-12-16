clear all
syms x1 x2 x3 x4
%% Setup
eps = 10^-8;
ops = sdpsettings('solver', 'mosek', 'debug', 1, 'verbose', 0);

Taylor4 = taylor(paperDyn([x1; x2; x3; x4]), [x1; x2; x3; x4], [0; 0; 0; 0], 'order', 4);
T4 = matlabFunction(Taylor4);
jac = jacobian(paperDyn([x1;x2;x3;x4]),[x1;x2;x3;x4]);
A = double(subs(jac,{x1,x2,x3,x4},{0,0,0,0}));
P = lyap(A',eye(length(A)));

sdpvar x_1 x_2 x_3 x_4
x = [x_1;x_2;x_3;x_4];
dyn = T4(x_1,x_2,x_3,x_4);
V = x'*P*x;
dV = jacobian(V,x);
sdpvar e1 e2
L1 = e1*x'*x;
L2 = e2*x'*x;
[s1,coefs_s1,~] = polynomial(x,2,0);
[s2,coefs_s2,~] = polynomial(x,4,0);
p = x'*x;

%%

gammaMin = .001;
gammaMax = 10;

for i = 1:20
    gamma = (gammaMin + gammaMax)/2;

    F = [sos(s2),
         sos(-(dV*dyn + L2 + s2*(gamma-V))),
         e2 >= eps];
 
    [sol,u,Q] = solvesos(F,[],ops,[coefs_s2;e2]);
    if sol.problem == 0
        gammaMin = gamma;
    else
        gammaMax = gamma;
    end
end

%%

betaMin = .001;
betaMax = 10;

for i = 1:20
    beta = (betaMin + betaMax)/2;

    F = [sos(-((V-gamma)+s1*(beta-p))),
     sos(s1)];
 
    [sol,u,Q] = solvesos(F,[],ops,[coefs_s1]);

    if sol.problem == 0
        betaMin = beta;
    else
        betaMax = beta;
    end
end

%%

[V,coefs_V,~] = polynomial(x,2,0);
dV = jacobian(V,x);

F = [sos(V - L1),
     sos(-(dV*dyn + L2 + s2*(gamma-V))),
     sos(-((V-gamma) +s1*(beta-p))),
     e1 >= eps];
     
[sol,u,Q] = solvesos(F,[],ops,[coefs_V; e1]);     