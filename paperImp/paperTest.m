clear all
beta = .2;
sdpvar z1 z2 z3 z4 z5 z6
sdpvar e11 e12 e13 e14 e15 e16
sdpvar e21 e22 e23 e24 e25 e26
z = [z1;z2;z3;z4;z5;z6];
fz = [z3 - z2*z3;
      z1*z3;
      0.4996*z4 - 0.4*z3 - 1.4994*z1 - .02*z5 + 0.02*z1*z4 + 0.4996*z1*z5 - 0.4996*z2*z4 + .02*z2*z5;
      z6 - z5*z6;
      z4*z6;
      0.4996*z1 + .02*z2 - .9986*z4 + .05*z5 - .5*z6 - .02*z1*z4 - 0.4996*z1*z5 + .4996*z2*z4 - .02*z2*z5];
G1 = z1^2 + z2^2 - 2*z2;
G2 = z4^2 + z5^2 - 2*z5;
p = z1^2 + z2^2 + 2*z3^2 + z4^2 + z5^2 + 2*z6^2;
%L_1 = e11*z1^2 + e12*z2^2 + e13*z3^2 + e14*z4^2 + e15*z5^2 + e16*z6^2;
%L_2 = e21*z1^2 + e22*z2^2 + e23*z3^2 + e24*z4^2 + e25*z5^2 + e26*z6^2;

[Vz,coefsV,~] = polynomial(z,4,2); %Define Barrier function to be a 6th order polynomial
[s2,coefs_s2,~] = polynomial(z,0,0); % s2
[s6,coefs_s6,~] = polynomial(z,1,0); % s6
[L_1,coefs_L_1,~] = polynomial(z,2,0); % L1
[L_2,coefs_L_2,~] = polynomial(z,2,0); % L2
[lam_11,coefs_Lam11,~] = polynomial(z,0,0); % lamb11
[lam_12,coefs_Lam12,~] = polynomial(z,0,0); % lamb11
[lam_21,coefs_Lam21,~] = polynomial(z,1,0); % lamb21
[lam_22,coefs_Lam22,~] = polynomial(z,1,0); % lamb21

ops = sdpsettings('solver', 'mosek');

dVz = jacobian(Vz,z); %derivatives
eps = 10^-8;

F = [sos((-s2*(beta-p))+Vz-(lam_11*G1)-(lam_12*G2)-L_1), 
     sos((-s6*(beta-p))-(dVz*fz)-(lam_21*G1)-(lam_22*G2)-L_2),
     sos(s2),sos(s6),
     %e11 >= eps, e12 >= eps, e13 >= eps, e14 >= eps, e15 >= eps, e16 >= eps,
     %e21 >= eps, e22 >= eps, e23 >= eps, e24 >= eps, e25 >= eps, e26 >= eps,
     %e11 + e12 + e13 + e14 + e15 + e16 >= .01,
     %e21 + e22 + e23 + e24 + e25 + e26 >= .01
     sos(L_1), sos(L_2)
     ];

 
solvesos(F,[],ops,[coefsV;coefs_s2;coefs_s6;coefs_Lam11;coefs_Lam21;coefs_Lam12;coefs_Lam22;
         %e11;e12;e13;e14;e15;e16;e21;e22;e23;e24;e25;e26
         coefs_L_1;coefs_L_2
         ])