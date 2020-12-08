function [dyn] = sym_gen_2()
%% function to generate the symbolic expressions required by SQP algorithm using CasADi
% addpath('E:\Program Files\MATLAB\CasADi')
addpath('C:\Users\pwest\Downloads\DCC - cleaned code\DCC - cleaned code\casadi')
import casadi.*;
%% Create stage functions
xk = SX.sym('xk',6);
par = SX.sym('par',12);
% create fuel cell dynamics
xdotk = grid_model(1, xk, par);

A = jacobian(xdotk,xk);
% create casadi function for the discrete dynamicsA = Function('f_x',{xk, uk, par}, {A});
A = Function('f_x',{xk, par}, {A});
dyn.A = A;



end
