function p_cv = p_lh_cv_interband(lh_energy, c_energy, Ep_1, Ep_2, Mb, Eg, m_eff_c, m_eff_lh)

mu = (abs(1/m_eff_c) + abs(1/m_eff_lh))^-1;
denominator = mu/lh_energy.*(Ep_1 + Ep_2 - abs(c_energy) - abs(lh_energy) - Eg)/m_eff_lh;

p_cv_squared = Mb^2*1/2.*(1 + 3./(1 + denominator));
p_cv_squared(denominator <= 0) = 0;
p_cv = sqrt(p_cv_squared);


end

