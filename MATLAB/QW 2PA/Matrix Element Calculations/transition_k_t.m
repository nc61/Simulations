function [k, step_function] = transition_k_t(E_pump, E_probe, Ec_final, Eh_initial, mu, Eg)

hbar = 6.626e-34/(2*pi);
step_function = ones(size(E_probe));

k_squared = 2*mu/hbar^2.*(E_pump + E_probe - Eg - (Ec_final + Eh_initial));
step_function(k_squared < 0) = 0;
k_squared(k_squared < 0) = 0;
k = sqrt(k_squared);


end

