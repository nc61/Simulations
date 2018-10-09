function [t, rho] = pulse_block(start, duration, R0, delta, prev_block, lambda_mat, gamma_mat)
% Duration - duration of the pulse
% R0

t_range = [start start+duration];

R0t = t_range;
delta_mat = delta;
ic = prev_block;


[t,rho] = ode45(@(t,rho)density_eqn_motion(t,rho,R0t,R0,delta_mat, gamma_mat, lambda_mat), t_range, ic);


end

