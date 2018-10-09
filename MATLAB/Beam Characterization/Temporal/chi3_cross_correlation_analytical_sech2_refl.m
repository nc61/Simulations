function normalized_transmission = chi3_cross_correlation_analytical_sech2_refl(delay_fs, delay_shift_fs, pump_refractive_index, pump_group_index, probe_group_index, pump_spot_size_m_HWeM, ...
    probe_spot_size_m_HWeM, pump_pulsewidth_fs_HWeM, probe_pulsewidth_fs_HWeM, sample_thickness_m, pump_energy_J, alpha_2_cm_per_GW, reflection_factor, amplitude_offset)

pump_reflectivity = ((pump_refractive_index - 1)/(pump_refractive_index + 1))^2;
pump_energy_J = (1 - pump_reflectivity)*pump_energy_J;
alpha_2_m_per_W = alpha_2_cm_per_GW*1e-11;
tau_gvd_s = sample_thickness_m/3e8*(pump_group_index - probe_group_index);
tau_gvd_s_reflection = sample_thickness_m/3e8*(-pump_group_index - probe_group_index);
tau_pump_s = 1e-15*pump_pulsewidth_fs_HWeM;
tau_probe_s = 1e-15*probe_pulsewidth_fs_HWeM;

reflected_energy = pump_reflectivity*pump_energy_J;

delay_fs = delay_fs - delay_shift_fs;
delay_s = 1e-15*delay_fs;
delay_s_reflection = delay_s + tau_gvd_s;
scale = 5;
t = linspace(-scale*tau_pump_s, scale*tau_pump_s, 200);

if (pump_group_index == probe_group_index)
    error('have not implemented the case for no GVD');
else
    E0 = 2*alpha_2_m_per_W*sample_thickness_m*pump_energy_J/(2*tau_pump_s)*1/(2*pi)*1/(pump_spot_size_m_HWeM^2 + probe_spot_size_m_HWeM^2)*1/tau_gvd_s;
    E0_reflection = reflection_factor*2*alpha_2_m_per_W*sample_thickness_m*reflected_energy/(2*tau_pump_s)*1/(2*pi)*1/(pump_spot_size_m_HWeM^2 + probe_spot_size_m_HWeM^2)*1/tau_gvd_s_reflection;

    cross_signal = zeros(size(delay_s));
    reflection = zeros(size(delay_s));

    for ind = 1:length(delay_fs)
        int = E0*sech(t/tau_pump_s).^2.*(tanh((t - delay_s(ind) - tau_gvd_s)/tau_probe_s) - tanh((t - delay_s(ind))/tau_probe_s));
        %int = 2/sqrt(pi)*E0*exp(-t.^2./tau_pump_s^2).*(tanh((t - delay_s(ind) - tau_gvd_s)/tau_probe_s) - tanh((t - delay_s(ind))/tau_probe_s));
        int_reflection = E0_reflection*sech(t/tau_pump_s).^2.*(tanh((t - delay_s_reflection(ind) - tau_gvd_s_reflection)/tau_probe_s) - tanh((t - delay_s_reflection(ind))/tau_probe_s));
        %int_reflection = 2/sqrt(pi)*E0_reflection*exp(-t.^2./tau_pump_s^2).*(tanh((t - delay_s_reflection(ind) - tau_gvd_s_reflection)/tau_probe_s) - tanh((t - delay_s_reflection(ind))/tau_probe_s));
        cross_signal(ind) = trapz(t, int);
        reflection(ind) = trapz(t, int_reflection);
    end
    normalized_transmission = cross_signal + reflection;
end

normalized_transmission = normalized_transmission + amplitude_offset;
