function [dyn] = Lyap_sym_gen()
%% function to generate the symbolic expressions required by SQP algorithm using CasADi
% addpath('E:\Program Files\MATLAB\CasADi')
addpath('C:\Users\pwest\Downloads\DCC - cleaned code\DCC - cleaned code\casadi')
import casadi.*;
%% Create stage functions
xk = SX.sym('xk',2);
Powers = SX.sym('p',3);
SEP = SX.sym('SEP',6);
% create fuel cell dynamics
[Lyapunov,~,~] = lyapunov_2(xk,SEP,Powers);

% A = jacobian(xdotk,xk);
% % create casadi function for the discrete dynamicsA = Function('f_x',{xk, uk, par}, {A});
% A = Function('f_x',{xk, par}, {A});
% dyn.A = A;



end
