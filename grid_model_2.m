function eqns = grid_model_2(t,x,u,p)
k = (-1)^(1/2);
B12 = -imag(1/(1.2*k));
B21 = B12;
B13 = -imag(1/(1*k));
B31 = B13;
B23 =  -imag(1/(0.8*k));
B32 = B23;
P1 = u(1); M1 = p(1); D1 = p(4);
P2 = u(2); M2 = p(2); D2 = p(5);
P3 = u(3); M3 = p(3); D3 = p(6);
V1 = u(4);
V2 = u(5);
V3 = u(6);
theta_12 = x(1);
theta_13 = x(2);
theta_23 = x(3);
w_1 = x(4);
w_2 = x(5);
w_3 = x(6);

theta_12_dot = w_1 - w_2;
theta_13_dot = w_1 - w_3;
theta_23_dot = w_2 - w_3;
w_1_dot = (1/M1)*(P1-((V1*V2*B12*sin(theta_12))+(V1*V3*B13*sin(theta_13)))-D1*w_1);
w_2_dot = (1/M2)*(P2-((V2*V1*B21*sin(-theta_12))+(V2*V3*B23*sin(theta_23)))-D2*w_2);
w_3_dot = (1/M3)*(P3-((V3*V1*B31*sin(-theta_13))+(V3*V2*B32*sin(-theta_23)))-D3*w_3);

eqns = [theta_12_dot;theta_13_dot;theta_23_dot;w_1_dot;w_2_dot;w_3_dot];

end

