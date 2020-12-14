% SOS formulation
%% Question 2 part 2
% SOS program for 2nd set of intial conditons
% Question 2, Homework 6
%Define parameter values
beta = 0.2;
Pm = 1.5;
D = 0.5;
P_max = 2.0;
SEP = 0.8481;
ops = sdpsettings('solver', 'sedumi');
sdpvar z1 z2 z3
sdpvar e1 e2 e3 e4 e5 e6
z = [z1;z2;z3];
fz = [(z1 - z2*z3);z1*z3;(Pm-P_max((z1*cos(SEP))+(1-z2)*sin(SEP))-D*z3)];
p = z1^2 + z2^2 + 2.0*z3^2;
L_1 = e1*z1^2 + e2*z2^2 + e3*z3^2;
L_2 = e4*z1^2 + e5*z2^2 + e6*z3^2;
G1 = z1^2 + z2^2 - 2*z2;

[Vz,coefsV,~] = polynomial(z,4,4); %Define Barrier function to be a 6th order polynomial

[s2,coefs_s2,~] = polynomial(z,2,2); % s2

[s6,coefs_s6,~] = polynomial(z,2,2); % s6

[lambda11,coefs_L1,~] = polynomial(z,0,0); % lamb11

[lambda21,coefsd_L2,~] = polynomial(z,1,1); % lamb21

dVz = jacobian(Vz,z); %derivatives


F = [sos((-s2*(beta-p))+Vz-(lambda11*G1)-L1), sos((-s6*(beta-p))-(dVz*fz)-(lambda21*G1)-L2),e1+e2+e3>0.01,e4+e5+e6>0.01,...
      e1>0,e2>0,e3>0,e4>0,e5>0,e6>0,sos(s2),sos(s6)]; 

solvesos(F,[],ops,[coefsV;coefs_s2;coefs_s6;coefs_L1;coefsd_L2;e1;e2;e3;e4;e5;e6])