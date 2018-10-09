hbar = 6.626e-34/(2*pi);
a = 2.83e-10;
kx = linspace(-pi/a, pi/a,100);
ky = linspace(-pi/a, pi/a,100);
kz = linspace(-pi/a, pi/a,100);
Eg_eV = 1.424;
E_opt_eV = 0.8;
gamma_eV = 0.2;
Eg_J = 1.6e-19*Eg_eV;
gamma_J = 1.6e-19*gamma_eV;
E_opt_J = 1.6e-19*E_opt_eV;

A0 = 1;
m_0 = 9.109e-31;
m_c = 0.067*m_0;
m_v = 1.2*m_0;

[Kx, Ky, Kz] = ndgrid(kx, ky, kz);

integrand = (1/m_0)^2*1.6e-19^2*1/((Eg_J + hbar^2*(Kx.^2 + Ky.^2 + Kz.^2)*(1/m_c + 1/m_v)).*(Eg_J + hbar^2*(Kx.^2 + Ky.^2 + Kz.^2)*(1/m_c + 1/m_v) - E_opt_J + 1i*gamma_J));
x_int =  squeeze(trapz(kz, trapz(ky, trapz(kx, integrand, 1), 2), 3));



