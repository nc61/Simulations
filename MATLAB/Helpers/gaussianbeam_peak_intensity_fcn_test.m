energy_J = 15.5e-6;
pulsewidth_fs_HW_inve_max = 100;
waist_mm_HW_inve_max = 1.2;

peak_intensity_GW_per_cm2 = gaussianbeam_peak_intensity_fcn(energy_J, waist_mm_HW_inve_max, pulsewidth_fs_HW_inve_max);

x_mm = linspace(-3*waist_mm_HW_inve_max, 3*waist_mm_HW_inve_max, 50);
y_mm  = linspace(-3*waist_mm_HW_inve_max, 3*waist_mm_HW_inve_max, 50);

dx_mm = x_mm(2) - x_mm(1);
t_fs = linspace(-3*pulsewidth_fs_HW_inve_max, 3*pulsewidth_fs_HW_inve_max, 50);
dt_fs = t_fs(2) - t_fs(1);

[X,Y,T] = ndgrid(x_mm,y_mm,t_fs);
phi = sqrt(peak_intensity_GW_per_cm2)*exp(-(X.^2 + Y.^2)/(2*waist_mm_HW_inve_max^2)).*exp(-T.^2/(2*pulsewidth_fs_HW_inve_max^2));
in1 = squeeze(trapz(1e-1*x_mm, 1e9*abs(phi).^2, 1));
in2 = squeeze(trapz(1e-1*y_mm, in1, 1));
out = squeeze(trapz(1e-15*t_fs, in2))



