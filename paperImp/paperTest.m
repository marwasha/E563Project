%% Part 1
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
L_1 = e11*z1^2 + e12*z2^2 + e13*z3^2 + e14*z4^2 + e15*z5^2 + e16*z6^2;
L_2 = e21*z1^2 + e22*z2^2 + e23*z3^2 + e24*z4^2 + e25*z5^2 + e26*z6^2;

[Vz,coefsV,~] = polynomial(z,2,1); %Define 2nd order Lyapunov Function
[s2,coefs_s2,~] = polynomial(z,2,0); % s2
[s6,coefs_s6,~] = polynomial(z,2,0); % s6
L1_coefs = coefficients(L_1,z);
L2_coefs = coefficients(L_2,z);
[lam_11,coefs_Lam11,~] = polynomial(z,0,0); % lamb11
[lam_12,coefs_Lam12,~] = polynomial(z,0,0); % lamb11
[lam_21,coefs_Lam21,~] = polynomial(z,1,0); % lamb21
[lam_22,coefs_Lam22,~] = polynomial(z,1,0); % lamb21

ops = sdpsettings('solver', 'sedumi');

dVz = jacobian(Vz,z); %derivatives
eps = 10^-8;

F = [sos((-s2*(beta-p))+Vz-(lam_11*G1)-(lam_12*G2)-L_1), 
     sos((-s6*(beta-p))-(dVz*fz)-(lam_21*G1)-(lam_22*G2)-L_2),
     sos(s2),sos(s6),
     sum(L1_coefs)>=0.1,  sum(L2_coefs)>=0.1, L1_coefs>=0,L2_coefs>=0,
     
     ];

 
solvesos(F,[],ops,[coefsV;coefs_s2;coefs_s6;coefs_Lam11;coefs_Lam21;coefs_Lam12;coefs_Lam22;
         L1_coefs;L2_coefs
         ])

%% Part 2
% Region of Attraction Visualization (EQ.22 from the paper)

[s1,coefs_s1,~] = polynomial(z,2,0); % s1
[lam_11,coefs_Lam11,~] = polynomial(z,1,0); % lamb11
[lam_12,coefs_Lam12,~] = polynomial(z,1,0); % lamb12

c = 8;

F = [sos((-s1*(c-Vz))-(s2*(p-beta))-((c-Vz)*(p-beta))-(lam_11*G1)-(lam_12*G2)-((p-beta)^2)), 
     sos(s1)];
 
 solvesos(F,[],ops,[coefs_s1;coefs_Lam11;coefs_Lam12])
 
 % Lyapunov functions as a matlab function
 Lyap_func = @(Z) (Z(1)*value(coefsV(1)))+(Z(2)*value(coefsV(2)))+(Z(3)*value(coefsV(3)))...
                  +(Z(4)*value(coefsV(4)))+(Z(5)*value(coefsV(5))) +(Z(6)*value(coefsV(6)))...
                  +((Z(1)^2)*value(coefsV(7)))+ (Z(1)*Z(2)*value(coefsV(8))) + ((Z(2)^2)*value(coefsV(9)))...
                  + (Z(1)*Z(3)*value(coefsV(10))) + (Z(2)*Z(3)*value(coefsV(11)))...
                  + ((Z(3)^2)*value(coefsV(12))) + (Z(1)*Z(4)*value(coefsV(13)))...
                  + (Z(2)*Z(4)*value(coefsV(14))) + (Z(3)*Z(4)*value(coefsV(15)))...
                  + ((Z(4)^2)*value(coefsV(16))) + (Z(1)*Z(5)*value(coefsV(17)))...
                  + (Z(2)*Z(5)*value(coefsV(18))) + (Z(3)*Z(5)*value(coefsV(19)))...
                  + (Z(4)*Z(5)*value(coefsV(20))) + ((Z(5)^2)*value(coefsV(21)))...
                  + (Z(1)*Z(6)*value(coefsV(22))) + (Z(2)*Z(6)*value(coefsV(23)))...
                  + (Z(3)*Z(6)*value(coefsV(24))) + (Z(4)*Z(6)*value(coefsV(25)))...
                  + (Z(5)*Z(6)*value(coefsV(26))) + ((Z(6)^2)*value(coefsV(27)));
                  
Lyap_input = @(X) [sin(X(1));1-cos(X(1));X(2);sin(X(3));1-cos(X(3));X(4)];
Test_point = [0,0,0,0];
Lyap_func(Lyap_input(Test_point))

 
 
 
 