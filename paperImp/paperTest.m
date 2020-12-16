%% Initalize V and beta
clear all
beta = .4;
eps = 10^-8;
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

ops = sdpsettings('solver', 'sedumi', 'verbose', 0);

dVz = jacobian(Vz,z); %derivatives
eps = 10^-8;

F = [sos((-s2*(beta-p))+Vz-(lam_11*G1)-(lam_12*G2)-L_1), 
     sos((-s6*(beta-p))-(dVz*fz)-(lam_21*G1)-(lam_22*G2)-L_2),
     sos(s2),sos(s6),
     sum(L1_coefs)>=0.1,  sum(L2_coefs)>=0.1, L1_coefs>=0,L2_coefs>=0,
     
     ];

 
solvesos(F,[],ops,[coefsV;coefs_s2;coefs_s6;coefs_Lam11;coefs_Lam21;
                   coefs_Lam12;coefs_Lam22;L1_coefs;L2_coefs])

%% Common

[s6,coefs_s6,~] = polynomial(z,2,0); % s2
[s8,coefs_s8,~] = polynomial(z,2,0); % s6
[s9,coefs_s9,~] = polynomial(z,2,0); % s9
[lam_11,coefs_Lam11,~] = polynomial(z,1,0); % lamb11
[lam_12,coefs_Lam12,~] = polynomial(z,1,0); % lamb11
[lam_21,coefs_Lam21,~] = polynomial(z,1,0); % lamb11
[lam_22,coefs_Lam22,~] = polynomial(z,1,0); % lamb11
[lam_31,coefs_Lam31,~] = polynomial(z,1,0); % lamb21
[lam_32,coefs_Lam32,~] = polynomial(z,1,0); % lamb21
     
%% 1a MAX c

c_min = eps;
c_max = 10;  

i = 0;
imax = 12;

while i < imax
    c = (c_min+c_max)/2;
    
    F = [sos(s6),sos(s8),sos(s9),
         sos(-s6*(beta-p) - (lam_21*G1)-(lam_22*G2) - (Vz-c)),
         sos(-s8*(c-Vz) - s9*(dVz*fz) - (lam_31*G1)-(lam_32*G2) - L_2)];
     
     [sol,u,Q] = solvesos(F,[],ops,[coefs_s6;coefs_s8;coefs_s9;coefs_Lam21;coefs_Lam22;coefs_Lam31;coefs_Lam32]);
     
     if sol.problem == 0
         c_min = c;
     else
         c_max = c;
         if i == imax - 1
             i = i - 1;
         end
     end
     i = i+1;
end
 
%% 1b MAX beta

beta_max = 10;
beta_min = eps;

j = 0; 
jmax = 12;

while j < jmax
   beta = (beta_max+beta_min)/2; 
    
   F = [sos(s6),sos(s8),sos(s9),
        sos(-s6*(beta-p) - (lam_21*G1)-(lam_22*G2) - (Vz-c)),
        sos(-s8*(c-Vz) - s9*(dVz*fz) - (lam_31*G1)-(lam_32*G2) - L_2)];
    
    [sol,u,Q] = solvesos(F,[],ops,[coefs_s6;coefs_s8;coefs_s9;coefs_Lam21;coefs_Lam22;coefs_Lam31;coefs_Lam32]);
     
     if sol.problem == 0
         beta_min = beta;
     else
         beta_max = beta;
         if j == jmax - 1
             j = j - 1;
         end
     end
     j = j+1;
   
end

%% 2a MIN c - STUCK

F = [sos(Vz - (lam_11*G1)-(lam_12*G2) - L_1);
     sos(-s6*(beta-p) - (lam_21*G1)-(lam_22*G2) - (Vz-c));
     sos(-s8*(c-Vz) - s9*(dVz*fz) - (lam_31*G1)-(lam_32*G2) - L_2);
     sos(s6)];
 
 [sol,u,Q] = solvesos(F,[],ops,[coefsV;coefs_s6;coefs_Lam11;coefs_Lam12;coefs_Lam21;coefs_Lam22;coefs_Lam31;coefs_Lam32]);