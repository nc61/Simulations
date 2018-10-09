function [alpha_ND_2, Ep_pump_J, Ep_2, absorption_edge] = ND2PA_TM_finite_parabolic_2(wavelength_pump_um, well_thickness_m, energies_c_J, wavefunctions_c, energies_lh_J, wavefunctions_lh, energies_hh_J, wavefunctions_hh)


% Universal constants
eps0 = 8.85e-12;
c = 3e8;
m0 = 9.109e-31;
q = 1.6e-19;
T = 295;
Eg = q*(1.519 - 5.405e-4*T^2/(T + 204));
E_kane = 25.7*q;
hbar = 6.626e-34/(2*pi);
urbach_tail_ev = 0;

% Effective masses for z (perpendicular) and transverse directions (z is
% growth direction)
gamma_1 = 6.8;
gamma_2 = 1.9;
m_lh_z = m0*1/(gamma_1 + 2*gamma_2);
m_lh_t = m0*1/(gamma_1 - gamma_2);

m_hh_z = m0*1/(gamma_1 - 2*gamma_2);
m_hh_t = m0*1/(gamma_1 + gamma_2);

m_e_z = 0.0665*m0;
m_e_t = m_e_z;
% QW properties

% Define photon energies
npoints = 1000;
Ep_pump_J = 1.24/wavelength_pump_um*q*ones(1, npoints);
                                                                 
% First do the LH1 - C2 transition
transition_energy_LH1C2 = Eg + energies_c_J(2) + energies_lh_J(1);
transition_energy_LH2C1 = Eg + energies_c_J(1) + energies_lh_J(2);
transition_energy_HH2C1 = Eg + energies_c_J(1) + energies_hh_J(2);
transition_energy_HH1C2 = Eg + energies_c_J(2) + energies_hh_J(1);

if (length(energies_hh_J) > 2)
    transition_energy_HH3C1 = Eg + energies_c_J(1) + energies_hh_J(3);
end

absorption_edge = Eg + energies_c_J(1) + energies_lh_J(1);
E_range = linspace(absorption_edge, absorption_edge + Ep_pump_J(1) - urbach_tail_ev*q, npoints);
Ep_2 = E_range - Ep_pump_J;

% Transition matrix elements
Px = sqrt(E_kane*m0/2);
Mb = sqrt((1/3)*Px^2);

p_c_intraband = @(n_f, n_i)-1i*hbar*m0/m_e_z*sum(wavefunctions_c(n_f,:).*[0 diff(wavefunctions_c(n_i,:))]);
p_lh_intraband = @(n_f, n_i)-1i*hbar*m0/m_lh_z*sum(wavefunctions_lh(n_f,:).*[0 diff(wavefunctions_lh(n_i,:))]);
p_hh_intraband = @(n_f, n_i)-1i*hbar*m0/m_hh_z*sum(wavefunctions_hh(n_f,:).*[0 diff(wavefunctions_hh(n_i,:))]);

I = @(Ep)1/2.*index2(1.24./(Ep/q), 0).*c.*eps0.*(Ep/hbar).^2;

% Joint density of states
mu_e_lh_t = (1/m_e_t + 1/m_lh_t)^-1;
g_2D_lh = mu_e_lh_t./(pi*hbar^2*well_thickness_m);

mu_e_hh_t = (1/m_e_t + 1/m_hh_t)^-1;
g_2D_hh = mu_e_hh_t./(pi*hbar^2*well_thickness_m);
Gamma_damping = 0;

W_tot = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create equation for the four possible transitions from LH1 - C2
transition_1_LH1C2 = p_lh_cv_interband(energies_lh_J(1), energies_c_J(2), Ep_pump_J, Ep_2, Mb, Eg, m_e_t, m_lh_t).*p_c_intraband(2, 1)./(Ep_2 - (energies_c_J(2) - energies_c_J(1)) + 1i*Gamma_damping);
transition_2_LH1C2 = p_lh_cv_interband(energies_lh_J(1), energies_c_J(2), Ep_pump_J, Ep_2, Mb, Eg, m_e_t, m_lh_t).*p_c_intraband(2, 1)./(Ep_pump_J - (energies_c_J(2) - energies_c_J(1)) + 1i*Gamma_damping);
transition_3_LH1C2 = p_lh_cv_interband(energies_lh_J(1), energies_c_J(2), Ep_pump_J, Ep_2, Mb, Eg, m_e_t, m_lh_t).*p_lh_intraband(2, 1)./((energies_lh_J(2) - energies_lh_J(1)) + Ep_pump_J + 1i*Gamma_damping);
transition_4_LH1C2 = p_lh_cv_interband(energies_lh_J(1), energies_c_J(2), Ep_pump_J, Ep_2, Mb, Eg, m_e_t, m_lh_t).*p_lh_intraband(2, 1)./((energies_lh_J(2) - energies_lh_J(1)) + Ep_2 + 1i*Gamma_damping);

H_sq_tot_LH1C2 = q^4/(16*m0^4).*abs(transition_1_LH1C2 + transition_2_LH1C2 + transition_3_LH1C2 + transition_4_LH1C2).^2;
W_ND_LH1C2 = @(Ep_total)2*pi/hbar*H_sq_tot_LH1C2*g_2D_lh.*heaviside(Ep_total - transition_energy_LH1C2);

W_tot = W_tot + W_ND_LH1C2(E_range);
figure(6)
plot(1.24./(E_range/q), W_ND_LH1C2(E_range).*Ep_2./(2*I(Ep_pump_J).*I(Ep_2))*1e11, 'r'), hold on



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create equation for the four possible transitions from LH2 - C1
transition_1_LH2C1 = p_lh_cv_interband(energies_lh_J(2), energies_c_J(1), Ep_pump_J, Ep_2, Mb, Eg, m_e_t, m_lh_t).*p_lh_intraband(1, 2)./((energies_lh_J(2) - energies_lh_J(1)) - Ep_pump_J + 1i*Gamma_damping);
transition_2_LH2C1 = p_lh_cv_interband(energies_lh_J(2), energies_c_J(1), Ep_pump_J, Ep_2, Mb, Eg, m_e_t, m_lh_t).*p_lh_intraband(1, 2)./((energies_lh_J(2) - energies_lh_J(1)) - Ep_2 + 1i*Gamma_damping);
transition_3_LH2C1 = p_lh_cv_interband(energies_lh_J(2), energies_c_J(1), Ep_pump_J, Ep_2, Mb, Eg, m_e_t, m_lh_t).*p_c_intraband(1, 2)./((energies_c_J(1) - energies_c_J(2)) - Ep_2 + 1i*Gamma_damping);
transition_4_LH2C1 = p_lh_cv_interband(energies_lh_J(2), energies_c_J(1), Ep_pump_J, Ep_2, Mb, Eg, m_e_t, m_lh_t).*p_c_intraband(1, 2)./((energies_c_J(1) - energies_c_J(2)) - Ep_pump_J + 1i*Gamma_damping);

H_sq_tot_LH2C1 = q^4/(16*m0^4).*abs(transition_1_LH2C1 + transition_2_LH2C1 + transition_3_LH2C1 + transition_4_LH2C1).^2;
W_ND_LH2C1 = @(Ep_total)2*pi/hbar*H_sq_tot_LH2C1*g_2D_lh.*heaviside(Ep_total - transition_energy_LH2C1);
W_tot = W_tot + W_ND_LH2C1(E_range);

plot(1.24./(E_range/q), W_ND_LH2C1(E_range).*Ep_2./(2*I(Ep_pump_J).*I(Ep_2))*1e11, 'g'), hold on



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create equation for the four possible transitions from HH2 - C1
transition_1_HH2C1 = p_hh_cv_interband(energies_hh_J(2), energies_c_J(1), Ep_pump_J, Ep_2, Mb, Eg, m_e_t, m_hh_t).*p_hh_intraband(1, 2)./((energies_hh_J(2) - energies_hh_J(1)) - Ep_pump_J + 1i*Gamma_damping);
transition_2_HH2C1 = p_hh_cv_interband(energies_hh_J(2), energies_c_J(1), Ep_pump_J, Ep_2, Mb, Eg, m_e_t, m_hh_t).*p_hh_intraband(1, 2)./((energies_hh_J(2) - energies_hh_J(1)) - Ep_2 + 1i*Gamma_damping);
transition_3_HH2C1 = p_hh_cv_interband(energies_hh_J(2), energies_c_J(1), Ep_pump_J, Ep_2, Mb, Eg, m_e_t, m_hh_t).*p_c_intraband(1, 2)./((energies_c_J(1) - energies_c_J(2)) - Ep_2 + 1i*Gamma_damping);
transition_4_HH2C1 = p_hh_cv_interband(energies_hh_J(2), energies_c_J(1), Ep_pump_J, Ep_2, Mb, Eg, m_e_t, m_hh_t).*p_c_intraband(1, 2)./((energies_c_J(1) - energies_c_J(2)) - Ep_pump_J + 1i*Gamma_damping);

H_sq_tot_HH2C1 = q^4/(16*m0^4).*abs(transition_1_HH2C1 + transition_2_HH2C1 + transition_3_HH2C1 + transition_4_HH2C1).^2;
W_ND_HH2C1 = @(Ep_total)2*pi/hbar*H_sq_tot_HH2C1*g_2D_hh.*heaviside(Ep_total - transition_energy_HH2C1);
W_tot = W_tot + W_ND_HH2C1(E_range);

plot(1.24./(E_range/q), W_ND_HH2C1(E_range).*Ep_2./(2*I(Ep_pump_J).*I(Ep_2))*1e11, 'b'), hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create equation for the four possible transitions from HH1 - C2
transition_1_HH1C2 = p_hh_cv_interband(energies_hh_J(1), energies_c_J(2), Ep_pump_J, Ep_2, Mb, Eg, m_e_t, m_hh_t).*p_c_intraband(2, 1)./(Ep_2 - (energies_c_J(2) - energies_c_J(1)) + 1i*Gamma_damping);
transition_2_HH1C2 = p_hh_cv_interband(energies_hh_J(1), energies_c_J(2), Ep_pump_J, Ep_2, Mb, Eg, m_e_t, m_hh_t).*p_c_intraband(2, 1)./(Ep_pump_J - (energies_c_J(2) - energies_c_J(1)) + 1i*Gamma_damping);
transition_3_HH1C2 = p_hh_cv_interband(energies_hh_J(1), energies_c_J(2), Ep_pump_J, Ep_2, Mb, Eg, m_e_t, m_hh_t).*p_hh_intraband(2, 1)./((energies_hh_J(2) - energies_hh_J(1)) + Ep_pump_J + 1i*Gamma_damping);
transition_4_HH1C2 = p_hh_cv_interband(energies_hh_J(1), energies_c_J(2), Ep_pump_J, Ep_2, Mb, Eg, m_e_t, m_hh_t).*p_hh_intraband(2, 1)./((energies_hh_J(2) - energies_hh_J(1)) + Ep_2 + 1i*Gamma_damping);

H_sq_tot_HH1C2 = q^4/(16*m0^4).*abs(transition_1_HH1C2 + transition_2_HH1C2 + transition_3_HH1C2 + transition_4_HH1C2).^2;
W_ND_HH1C2 = @(Ep_total)2*pi/hbar*H_sq_tot_HH1C2*g_2D_hh.*heaviside(Ep_total - transition_energy_HH1C2);
W_tot = W_tot + W_ND_HH1C2(E_range);

plot(1.24./(E_range/q), W_ND_HH1C2(E_range).*Ep_2./(2*I(Ep_pump_J).*I(Ep_2))*1e11, 'k'), hold on


if length(energies_hh_J) > 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create equation for the four possible transitions from HH3 - C1
transition_1_HH3C1 = p_hh_cv_interband(energies_hh_J(3), energies_c_J(1), Ep_pump_J, Ep_2, Mb, Eg, m_e_t, m_hh_t).*p_hh_intraband(1, 3)./((energies_hh_J(3) - energies_hh_J(1)) - Ep_pump_J + 1i*Gamma_damping);
transition_2_HH3C1 = p_hh_cv_interband(energies_hh_J(3), energies_c_J(1), Ep_pump_J, Ep_2, Mb, Eg, m_e_t, m_hh_t).*p_hh_intraband(1, 3)./((energies_hh_J(3) - energies_hh_J(1)) - Ep_2 + 1i*Gamma_damping);

H_sq_tot_HH3C1 = q^4/(16*m0^4).*abs(transition_1_HH3C1 + transition_2_HH3C1).^2;
W_ND_HH3C1 = @(Ep_total)2*pi/hbar*H_sq_tot_HH3C1*g_2D_hh.*heaviside(Ep_total - transition_energy_HH3C1);
W_tot = W_tot + W_ND_HH3C1(E_range);

plot(1.24./(E_range/q), W_ND_HH3C1(E_range).*Ep_2./(2*I(Ep_pump_J).*I(Ep_2))*1e11, 'm'), hold on

end
legend('LH1C2', 'LH2C1', 'HH2C1', 'HH1C2')
xlabel('\lambda (\mu m)'), ylabel('\beta_2')
hold off

alpha_ND_2 = (W_tot).*Ep_2./(2*I(Ep_pump_J).*I(Ep_2))*1e11;


end


