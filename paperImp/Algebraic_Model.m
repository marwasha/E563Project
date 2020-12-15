function [fz,gz,gx,g_compare] = Algebraic_Model(p,x,z,SEP)
x1 = x(1);
x2 = x(2);
Pm = p(1);
P_max = p(2);
z1 = z(1);
z2 = z(2);
z3 = z(3);
fz = [z3 - z2*z3;
      z1*z3;
      Pm - P_max*(cos(SEP) + (1-z2)*sin(SEP)) - D*z3];
gz = z1^2+z2^2-2*z3;
g_compare = [sin(x1)-z1;1-cos(x1)-z2;x2-z3];
gx = (sin(x1)^2) + ((1-cos(x1))^2) + 2*x2;
end