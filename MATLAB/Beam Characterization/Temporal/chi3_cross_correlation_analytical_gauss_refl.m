function [normalized_transmission_total, normalized_transmission, normalized_transmission_reflection] = chi3_cross_correlation_analytical_gauss_refl(delay_fs, delay_shift_fs, pump_refractive_index, pump_group_index, probe_group_index, pump_spot_size_m_HWeM, ...
    probe_spot_size_m_HWeM, pump_pulsewidth_fs_HWeM, probe_pulsewidth_fs_HWeM, sample_thickness_m, pump_energy_J, alpha_2_cm_per_GW, reflection_factor, amplitude_offset, include_reflections)

pump_reflectivity = ((pump_refractive_index - 1)/(pump_refractive_index + 1))^2;
pump_energy_J = (1 - pump_reflectivity)*pump_energy_J;
alpha_2_m_per_W = alpha_2_cm_per_GW*1e-11;
tau_gvd_s = -sample_thickness_m/3e8*(pump_group_index - probe_group_index);
tau_gvd_s_reflection = sample_thickness_m/3e8*(-pump_group_index - probe_group_index);
pump_pulsewidth_s_HWeM = 1e-15*pump_pulsewidth_fs_HWeM;
probe_pulsewidth_s_HWeM = 1e-15*probe_pulsewidth_fs_HWeM;

reflected_energy_J = pump_reflectivity*pump_energy_J;

delay_fs = delay_fs - delay_shift_fs;
delay_s = -1e-15*delay_fs;
delay_s_reflection = delay_s + tau_gvd_s;

cross_correlation_width_s_HWeM = sqrt(pump_pulsewidth_s_HWeM^2 + probe_pulsewidth_s_HWeM^2);

if (pump_group_index == probe_group_index)
    error('have not implemented the case for no GVD');
else
    normalized_transmission = -2*alpha_2_m_per_W*pi^(-3/2)*pump_energy_J*sample_thickness_m*1/(pump_spot_size_m_HWeM^2 + probe_spot_size_m_HWeM^2) ...
        *sqrt(pi)/2/tau_gvd_s*(erf(delay_s/cross_correlation_width_s_HWeM) - erf((delay_s - tau_gvd_s)/cross_correlation_width_s_HWeM));
    
    if include_reflections == 1
    normalized_transmission_reflection = -reflection_factor*2*alpha_2_m_per_W*pi^(-3/2)*reflected_energy_J*sample_thickness_m*1/(pump_spot_size_m_HWeM^2 + probe_spot_size_m_HWeM^2)*sqrt(1/(pump_pulsewidth_s_HWeM^2 + probe_pulsewidth_s_HWeM^2)) ...
        *sqrt(pi)/2*cross_correlation_width_s_HWeM/tau_gvd_s_reflection*(erf(delay_s_reflection/cross_correlation_width_s_HWeM) - erf((delay_s_reflection - tau_gvd_s_reflection)/cross_correlation_width_s_HWeM));
    else
        normalized_transmission_reflection = 0;
    end
    
end

normalized_transmission = normalized_transmission + amplitude_offset;
normalized_transmission_reflection = normalized_transmission_reflection + amplitude_offset;
normalized_transmission_total = normalized_transmission + normalized_transmission_reflection;