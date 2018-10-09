function gamma_minv_Winv_ij = gamma_nd(overlap_parameter, wavelength_um, A_eff_um2_i, A_eff_um2_j, n2_m2_per_W, beta_2pa_cm_per_GW)

A_eff_m2_i = A_eff_um2_i*1e-12;
A_eff_m2_j = A_eff_um2_j*1e-12;
beta_2pa_m_per_W = beta_2pa_cm_per_GW*1e-9*1e-2;

k_0 = 2*pi./(wavelength_um*1e-6);
gamma_minv_Winv_ij = overlap_parameter/sqrt(A_eff_m2_i*A_eff_m2_j)*(k_0*n2_m2_per_W + 1i*beta_2pa_m_per_W);

end

