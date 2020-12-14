function eqns = grid_model_CORA(x,u)

dx = grid_model_2(1,x(1:6), u, x(7:12));

eqns = [dx; zeros(6,1)];

end