% SOS formulation
%% Question 2 part 2
% SOS program for 2nd set of intial conditons
% Question 2, Homework 6
%Define parameter values
beta = 1;
Pm = 1.5;
D = 0.5;
P_max = 2.0;
SEP = asin(Pm/P_max);
ops = sdpsettings('solver', 'mosek');
sdpvar z1 z2 z3
z = [z1;z2;z3];
fz = [z3 - z2*z3;
      z1*z3;
      (Pm-P_max*((z1*cos(SEP))+(1-z2)*sin(SEP))-D*z3)];

p = z1^2 + z2^2 + 2.0*z3^2;
G1 = z1^2 + z2^2 - 2*z2;

[Vz,coefsV,~] = polynomial(z,4,4); %Define Barrier function to be a 6th order polynomial
[s2,coefs_s2,~] = polynomial(z,2,2); % s2
[s6,coefs_s6,~] = polynomial(z,2,2); % s6
[L_1,coefs_L_1,~] = polynomial(z,2,2); % L1
[L_2,coefs_L_2,~] = polynomial(z,2,2); % L2
[lambda11,coefs_Lam1,~] = polynomial(z,2,2); % lamb11
[lambda21,coefs_Lam2,~] = polynomial(z,2,2); % lamb21

dVz = jacobian(Vz,z); %derivatives

eps = 10^-8;


F = [sos((-s2*(beta-p))+Vz-(lambda11*G1)-L_1), 
     sos((-s6*(beta-p))-(dVz*fz)-(lambda21*G1)-L_2),
     sos(s2),sos(s6),sos(L_1),sos(L_2)]; 

solvesos(F,-sum(coefsV),ops,[coefsV;coefs_s2;coefs_s6;coefs_Lam1;coefs_Lam2;coefs_L_1;coefs_L_2])

%% ROA

ops = sdpsettings('solver', 'mosek', 'debug', 1);

[s1,coefs_s1,~] = polynomial(z,2,2);
[lam,coefs_lam,~] = polynomial(z,2,2);
c = .1;

F2 = [sos(-s1*(c-Vz) - s2*(p-beta) - (c-Vz)*(p-beta) - lam*G1 - (p-beta)^2);
      sos(s1)];
solvesos(F2,[],ops,[coefs_s1;coefs_lam])  
