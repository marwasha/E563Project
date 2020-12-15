function [dyn] = Algebraic_Model_sym_gen()
addpath('ToolBoxes\casadi-windows-matlabR2016a-v3.5.5')
import casadi.*;
%% Create stage functions
x = SX.sym('xk',2);
z = SX.sym('zk',3);
p = SX.sym('p',2);
[fz,gz,gx,gcompare] = Algebraic_Model(p,x,z,SEP);
dfdz = jacobian(fz,z);
dfdz = Function('df_dz',{z, p}, {dfdz});
dyn.dfdz = dfdz;
dgdz = jacobian(gz,z);
dgdz = Function('dg_dz',{z, p}, {dgdz});
dyn.dgdz = dgdz;
dgdx = jacobian(gx,x);
dgdx = Function('dg_dx',{z, p,x}, {dgdx});
dyn.dgdx = dgdx;
dgcompdx = jacobian(gcompare,x);
dgcompdx = Function('dgcomp_dx',{z, p,x}, {dgcompdx});
dyn.dgcompdx = dgcompdx;
dgcompdz = jacobian(gcompare,z);
dgcompdz = Function('dcomp_dz',{z, p,x}, {dgcompdz});
dyn.dgcompdz = dgcompdz;



end