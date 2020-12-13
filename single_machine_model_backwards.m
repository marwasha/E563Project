function [eqns] = single_machine_model_backwards(t,x,p)
P_m = p(1);
D = p(3);
P = p(2);
eqns = [-x(2);-1*(P_m - P*sin(x(1))-D*x(2))];
end

