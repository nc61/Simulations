function [wavefunction, energy] = GaAs_infinite_well(d, n, band, npoints)

m0 = 9.109e-31;

hbar = 6.626e-34/(2*pi);

% Effective masses for z (perpendicular) and transverse directions (z is
% growth direction)
gamma_1 = 6.8;
gamma_2 = 1.9;
m_lh_z = m0*1/(gamma_1 + 2*gamma_2);
m_hh_z = m0*1/(gamma_1 - 2*gamma_2);
m_e_z = 0.0665*m0;

if (strcmp(band,'lh'))
    m_eff = m_lh_z;
elseif strcmp(band,'hh')
    m_eff = m_hh_z;
elseif (strcmp(band,'c'))
    m_eff = m_e_z;
else
    m_eff = 0;
end

z = linspace(0, d, npoints);

wavefunction = sqrt(2/d)*sin(n*pi*z./d);
energy = n.^2*pi^2*hbar^2/(2*m_eff*d^2);

end

