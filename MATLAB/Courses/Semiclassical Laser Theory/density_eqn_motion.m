function rho_dot = density_eqn_motion(t,rho,R0t,R0,delta, gamma_mat, lambda_mat)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


rho_dot = zeros(length(rho),1);

lambda_a = lambda_mat(1);
lambda_b = lambda_mat(2);

gamma_aa = gamma_mat(1);
gamma_bb = gamma_mat(2);
gamma_ab = gamma_mat(3);


for ind = 1:3:3*length(delta)-2
    rho_dot(ind) = -gamma_aa*(rho(ind) - lambda_a) - 1i*R0/2.*(rho(ind+2)-conj(rho(ind+2)));
    rho_dot(ind+1) = -gamma_bb*(rho(ind+1) - lambda_b)  + 1i*R0/2.*(rho(ind+2) - conj(rho(ind+2)));
    rho_dot(ind+2) = -(1i*delta((ind+2)/3) + (gamma_aa + gamma_bb)/2 + gamma_ab)*rho(ind+2) - 1i*R0/2*(rho(ind) - rho(ind+1));
end

