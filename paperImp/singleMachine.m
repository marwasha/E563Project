% SOS formulation
%% Question 2 part 2
% SOS program for 2nd set of intial conditons
% Question 2, Homework 6
%Define parameter values
clear all
beta = .4;
Pm = 1.5;
D = 0.5;
P_max = 2.0;
SEP = asin(Pm/P_max);
ops = sdpsettings('solver', 'sedumi');
sdpvar z1 z2 z3
sdpvar e11 e12 e13
sdpvar e21 e22 e23
z = [z1;z2;z3];
fz = [z3 - z2*z3;
      z1*z3;
      (Pm-P_max*((z1*cos(SEP))+(1-z2)*sin(SEP))-D*z3)];

p = z1^2 + z2^2 + z3^2;
G = z1^2 + z2^2 - 2*z2;
L_1 = e11*z1^2 + e12*z2^2 + e13*z3^2;
L_2 = e21*z1^2 + e22*z2^2 + e23*z3^2; 

[Vz,coefsV,~] = polynomial(z,2,1); %Define Barrier function to be a 6th order polynomial
[s2,coefs_s2,~] = polynomial(z,2,0); % s2
[s6,coefs_s6,~] = polynomial(z,2,0); % s6
L1_coefs = coefficients(L_1,z);
L2_coefs = coefficients(L_2,z);
[lam1,coefs_Lam1,~] = polynomial(z,2,0); % lamb11
[lam2,coefs_Lam2,~] = polynomial(z,2,0); % lamb21
dVz = jacobian(Vz,z); %derivatives

eps = 10^-8;

F = [sos((-s2*(beta-p))+Vz-(lam1*G)-L_1), 
     sos((-s6*(beta-p))-(dVz*fz)-(lam2*G)-L_2),
     sos(s2),sos(s6),
     sum(L1_coefs)>=0.1,sum(L2_coefs)>=0.1,
     L1_coefs>=0,L2_coefs>=0
     ]; 

solvesos(F,-sum(coefsV),ops,[coefsV;coefs_s2;coefs_s6;coefs_Lam1;
                             coefs_Lam2;L1_coefs;L2_coefs])

%% ROA

ops = sdpsettings('solver', 'mosek', 'debug', 1);

[s1,coefs_s1,~] = polynomial(z,2,2);
[lam,coefs_lam,~] = polynomial(z,2,2);
c = .1;

F2 = [sos(-s1*(c-Vz) - s2*(p-beta) - (c-Vz)*(p-beta) - lam*G - (p-beta)^2);
      sos(s1)];
solvesos(F2,[],ops,[coefs_s1;coefs_lam])  
