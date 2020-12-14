x0_NL = [0.3358   -0.0000    0.3029    0.0000   -0.0460   -0.0000]';
p = [0.5,1,1.2,0.5,0.5,0.5]';
u = [0.4,0.4,-0.8,1,1,1]';
uF = [0.4,0.0,-0.8,1,1,1]';

dyn = sym_gen_2;
A1 = full(dyn.A(x0_NL,u,p));
B1 = full(dyn.B(x0_NL,u,p));
F1 = B1(:,2:6)*u(2:6);
B1 = B1(:,1);

A2 = full(dyn.A(x0_NL,uF,p));
B2 = full(dyn.B(x0_NL,uF,p));
F2 = B2(:,2:6)*uF(2:6);
B2 = B2(:,1);