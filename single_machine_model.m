function [eqns] = single_machine_model(t,x,p,SEP)
P_m = p(1);
D = p(3);
P = p(2);
eqns = [x(2);(P_m - P*sin(x(1)+SEP)-D*x(2))];
end

