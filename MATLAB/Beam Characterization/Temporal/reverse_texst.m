pump_energy_J = 4e-6;
alpha_2_cm_per_GW = 1.43;
sample_thickness_m = 80e-6;
pump_spot_size_m_HWeM = 407e-6;
probe_spot_size_m_HWeM = 75.7e-6;
probe_pulsewidth_fs_HWeM = 113;
pump_pulsewidth_fs_HWeM = 82;
pump_group_index = -3.4441;
probe_group_index = 3.7620;
pump_refractive_index = 3.4258;


%normalized_transmission = chi3_cross_correlation_analytical(delay_fs, delay_shift_s, pump_refractive_index, pump_group_index, probe_group_index, pump_spot_size_m_HWeM, ...
 %   probe_spot_size_m_HWeM, pump_pulsewidth_fs_HWeM, probe_pulsewidth_fs_HWeM, sample_thickness_m, pump_energy_J, alpha_2_cm_per_GW);

fit_normalized_transmission = @(fit_parameters, delay_fs) chi3_cross_correlation_analytical_sech2(delay_fs, fit_parameters(3), pump_refractive_index, pump_group_index, probe_group_index, pump_spot_size_m_HWeM, ...
    probe_spot_size_m_HWeM, fit_parameters(2), probe_pulsewidth_fs_HWeM, sample_thickness_m, pump_energy_J, fit_parameters(1));

[scan_delay_fs, scan_norm_trans] = read_scan('plot', 'off', 'normalizeToProbe', 'on', 'probeReading', 6.2, 'probeSensitivity', 20, 'scanSensitivity', 20, 'centerDataToPeak', 'on', 'micrometerUnits', 'in');
scan_delay_fs = -1*scan_delay_fs;
fit_params = lsqcurvefit(fit_normalized_transmission, [1.55, 80, 200], scan_delay_fs + 1e15*sample_thickness_m/3e8*(pump_group_index - probe_group_index), scan_norm_trans);

%fit_params = [1.55, 80e-15, 0];
fit_results = fit_normalized_transmission(fit_params, scan_delay_fs);
plot(scan_delay_fs, fit_results, 'b')
hold on
plot(scan_delay_fs, scan_norm_trans, 'ro')
hold off