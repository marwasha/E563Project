function eqns = singleMachineDyn(t,x,u,p)

delta = x(1);
omega = x(2);
P_m = u;
P_max = p(1);
M = p(2);
D = P(3);

eqns = zeros(2,1);

eqns(1) = x(2);
eqns(2) = (P_m - P_max*sin(delta) - D*omega)/M;

end