function normalized_transmission= cross_correlation_analytical(pump_refractive_index, pump_group_index, probe_group_index, pump_spot_size_m_HWeM_size_m_HWeM, ...
    probe_spot_size_m_HWeM, pump_pulsewidth_s_HWeM, probe_pulsewidth_s_HWeM, sample_thickness_m, pump_energy_J, alpha_2_cm_per_GW)

pump_reflectivity = 0.31;
pump_energy_J = (1 - pump_reflectivity)*3.10e-6;
alpha_2_cm_per_GW = 1.43;
alpha_2_m_per_W = alpha_2_cm_per_GW*1e-11;
sample_thickness_m = 280e-6;
pump_spot_size_m_HWeM = 407e-6;
probe_spot_size_m_HWeM = 75.7e-6;
probe_pulsewidth_s_HWeM = 138e-15;
pump_pulsewidth_s_HWeM = 82e-15;
pump_group_index = 3.4441;
probe_group_index = 3.7620;
tau_gvd = sample_thickness_m/3e8*(pump_group_index - probe_group_index);

cross_correlation_width_s_HWeM = sqrt(pump_pulsewidth_s_HWeM^2 + probe_pulsewidth_s_HWeM^2);

delay_scale = 8;
delay_s = delay_scale*linspace(-probe_pulsewidth_s_HWeM, probe_pulsewidth_s_HWeM, 500);

if (pump_group_index == probe_group_index)
    
    normalized_transmission = -2*alpha_2_m_per_W*pi^(-3/2)*pump_energy_J*sample_thickness_m*1/(pump_spot_size_m_HWeM^2 + probe_spot_size_m_HWeM^2)*sqrt(1/(pump_pulsewidth_s_HWeM^2 + probe_pulsewidth_s_HWeM^2)) ...
        *exp(-delay_s.^2/(pump_pulsewidth_s_HWeM^2 + probe_pulsewidth_s_HWeM^2));
    
else
    
    normalized_transmission_GVD = -2*alpha_2_m_per_W*pi^(-3/2)*pump_energy_J*sample_thickness_m*1/(pump_spot_size_m_HWeM^2 + probe_spot_size_m_HWeM^2)*sqrt(1/(pump_pulsewidth_s_HWeM^2 + probe_pulsewidth_s_HWeM^2)) ...
        *sqrt(pi)/2*cross_correlation_width_s_HWeM/tau_gvd*(erf(delay_s/cross_correlation_width_s_HWeM) - erf((delay_s - tau_gvd)/cross_correlation_width_s_HWeM));
end

figure(1)
plot(delay_s, normalized_transmission, 'r');
hold on
plot(delay_s, normalized_transmission_GVD, 'bo');
hold off
