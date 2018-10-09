wavelength_probe_um = 1.25;
filename_1 = 'TM1_1900.mat';
mode_numer_1 = 1;
polarization_1 = 'z';

wavelength_pump_um = 1.9;
filename_2 = 'TM1_1250.mat';
mode_number_2 = 1;
polarization_2 = 'z';

n2_m2_per_W = 0;
n2_m2_per_W_pump = 1e-15;
beta_2pa_cm_per_GW = 25;
beta_2pa_cm_per_GW_pump = 0;

[A_eff_um2_probe, A_eff_um2_pump, overlap] = mode_area_and_overlap(filename_1, mode_numer_1, polarization_1, filename_2, mode_number_2, polarization_2)
gamma_pe = gamma_nd(overlap, wavelength_probe_um, A_eff_um2_probe, A_eff_um2_pump, n2_m2_per_W, beta_2pa_cm_per_GW)
gamma_ee = gamma_nd(overlap, wavelength_probe_um, A_eff_um2_pump, A_eff_um2_pump, n2_m2_per_W_pump, beta_2pa_cm_per_GW_pump)

[n_eff_pump, beta_1_fs_per_mm_pump, beta_2_fs2_per_mm_pump] = neff_dispersion('QW_dispersion.mat', wavelength_pump_um)
[n_eff_probe, beta_1_fs_per_mm_probe, beta_2_fs2_per_mm_probe] = neff_dispersion('QW_dispersion.mat', wavelength_probe_um) 

gvm_parameter_fs_per_mm = beta_1_fs_per_mm_probe - beta_1_fs_per_mm_pump

