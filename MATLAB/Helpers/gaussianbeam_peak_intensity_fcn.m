function peak_intensity_GW_per_cm2 = gaussianbeam_peak_intensity_fcn(energy_J, waist_mm_HW_inve_max, pulsewidth_fs_HW_inve_max)

peak_intensity_GW_per_cm2 = (1e-9)*energy_J/(pi^(3/2)*(1e-1*waist_mm_HW_inve_max)^2*(1e-15)*pulsewidth_fs_HW_inve_max);

end

