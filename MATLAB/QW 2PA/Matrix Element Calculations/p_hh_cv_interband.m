function p_cv = p_hh_cv_interband(hh_energy, c_energy, Ep_1, Ep_2, Mb, Eg, m_eff_c, m_eff_hh)

mu = (abs(1/m_eff_c) + abs(1/m_eff_hh))^-1;
denominator = mu/hh_energy.*(Ep_1 + Ep_2 - abs(c_energy) - abs(hh_energy) - Eg)/m_eff_hh;

p_cv_squared = Mb^2*3/2.*(1 - 1./(1 + denominator));
p_cv_squared(denominator <= 0) = 0;
p_cv = sqrt(p_cv_squared);


end

