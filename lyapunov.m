function [Lyapunov,PE,KE] = lyapunov(x,SEP,Powers)
k = (-1)^(1/2);
B12 = -imag(1/(1.2*k));
B21 = B12;
B13 = -imag(1/(1*k));
B31 = B13;
B23 =  -imag(1/(0.8*k));
B32 = B23;
M1 = 0.5;
M2 = 1.0;
M3 = 1.2;
P1 = Powers(1);
P2 = Powers(2);
P3 = Powers(3);
Lyapunov =  1/2*[x(2) x(4) x(6)]*diag([M1 M2 M3])*[x(2);x(4);x(6)] - ... 
                1/2*(B12*(cos(x(1)-x(3))-cos(SEP(1)-SEP(3)))+B13*(cos(x(1)-x(5))-cos(SEP(1)-SEP(5))) + ... 
                B21*(cos(x(3)-x(1))-cos(SEP(3)-SEP(1)))+B23*(cos(x(3)-x(5))-cos(SEP(3)-SEP(5))) + ...
                B31*(cos(x(5)-x(1))-cos(SEP(5)-SEP(1)))+B32*(cos(x(5)-x(3))-cos(SEP(5)-SEP(3))))- ...
                (P1*(x(1)-SEP(1))+P2*(x(3)-SEP(3))+P3*(x(5)-SEP(5)));
            
PE =        -1/2*(B12*(cos(x(1)-x(3))-cos(SEP(1)-SEP(3)))+B13*(cos(x(1)-x(5))-cos(SEP(1)-SEP(5))) + ...
                B21*(cos(x(3)-x(1))-cos(SEP(3)-SEP(1)))+B23*(cos(x(3)-x(5))-cos(SEP(3)-SEP(5))) + ...
                B31*(cos(x(5)-x(1))-cos(SEP(5)-SEP(1)))+B32*(cos(x(5)-x(3))-cos(SEP(5)-SEP(3))))- ...
                (P1*(x(1)-SEP(1))+P2*(x(3)-SEP(3))+P3*(x(5)-SEP(5)));
            
KE =         1/2*[x(2) x(4) x(6)]*diag([M1 M2 M3])*[x(2);x(4);x(6)]; 
end

