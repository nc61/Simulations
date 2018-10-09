function normalized_transmission = nlr_cross_correlation_analytical(delay_fs, delay_shift_fs, pump_refractive_index, pump_group_index, probe_group_index, pump_spot_size_m_HWeM, ...
    probe_spot_size_m_HWeM, pump_pulsewidth_fs_HWeM, probe_pulsewidth_fs_HWeM, sample_thickness_m, pump_energy_J, n_2_m2_per_W, amplitude_offset)

pump_reflectivity = ((pump_refractive_index - 1)/(pump_refractive_index + 1))^2;
pump_energy_J = (1 - pump_reflectivity)*pump_energy_J;
tau_gvd_fs = 1e15*sample_thickness_m/3e8*(pump_group_index - probe_group_index);

cross_correlation_width_fs_HWeM = sqrt(pump_pulsewidth_fs_HWeM^2 + probe_pulsewidth_fs_HWeM^2);
delay_fs = delay_fs - delay_shift_fs;

if (pump_group_index == probe_group_index)
    normalized_transmission = -2i*k_0_minv*n_2_m2_per_W*pi^(-3/2)*pump_energy_J*sample_thickness_m*1/(pump_spot_size_m_HWeM^2 + probe_spot_size_m_HWeM^2)*1e15*sqrt(1/(pump_pulsewidth_fs_HWeM^2 + probe_pulsewidth_fs_HWeM^2)) ...
        *exp(-delay_fs.^2/(pump_pulsewidth_fs_HWeM^2 + probe_pulsewidth_fs_HWeM^2));
else
    normalized_transmission = -2i*k_0_minv*n_2_m2_per_W*pi^(-3/2)*pump_energy_J*sample_thickness_m*1/(pump_spot_size_m_HWeM^2 + probe_spot_size_m_HWeM^2)*1e15*sqrt(1/(pump_pulsewidth_fs_HWeM^2 + probe_pulsewidth_fs_HWeM^2)) ...
        *sqrt(pi)/2*cross_correlation_width_fs_HWeM/tau_gvd_fs*(erf(delay_fs/cross_correlation_width_fs_HWeM) - erf((delay_fs - tau_gvd_fs)/cross_correlation_width_fs_HWeM));
end

normalized_transmission = normalized_transmission + amplitude_offset;
