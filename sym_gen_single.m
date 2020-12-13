function [dyn] = sym_gen_single()
%% function to generate the symbolic expressions required by SQP algorithm using CasADi
% addpath('E:\Program Files\MATLAB\CasADi')
addpath('C:\Users\pwest\Downloads\DCC - cleaned code\DCC - cleaned code\casadi')
import casadi.*;
%% Create stage functions
xk = SX.sym('xk',2);
par = SX.sym('par',3);
% create fuel cell dynamics
xdotk = single_machine_model(1, xk, par);

A = jacobian(xdotk,xk);
% create casadi function for the discrete dynamicsA = Function('f_x',{xk, uk, par}, {A});
A = Function('f_x',{xk, par}, {A});
dyn.A = A;



end
