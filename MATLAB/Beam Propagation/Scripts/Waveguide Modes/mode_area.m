

load('mode_1300.mat', 'mode1_Ey');
load('mode_1300.mat', 'x', 'y', 'z');

mode_1300 = mode1_Ey;

load('mode_1800.mat', 'mode1_Ey');
mode_1800 = mode1_Ey;

dy = [y(2) - y(1); diff(y)];
dz = [z(2) - z(1); diff(z)];

mode_irradiance_1300 = abs(mode_1300).^2;
mode_irradiance_squared_1300 = abs(mode_1300).^4;
irradiance_integral_1300 = trapz(z, trapz(y, mode_irradiance_1300, 2), 3);
irradiance_squared_integral_1300 = trapz(z, trapz(y, mode_irradiance_squared_1300, 2), 3);
mode_area_m2_1300 = irradiance_integral_1300^2/irradiance_squared_integral_1300;
mode_area_um2_1300 = mode_area_m2_1300*1e12

mode_irradiance_1800 = abs(mode_1800).^2;
mode_irradiance_squared_1800 = abs(mode_1800).^4;
irradiance_integral_1800 = trapz(z, trapz(y, mode_irradiance_1800, 2), 3);
irradiance_squared_integral_1800 = trapz(z, trapz(y, mode_irradiance_squared_1800, 2), 3);
mode_area_m2_1800 = irradiance_integral_1800^2/irradiance_squared_integral_1800;
mode_area_um2_1800 = mode_area_m2_1800*1e12

overlap_integral = trapz(z, trapz(y, mode_irradiance_1300.*mode_irradiance_1800, 2), 3);
overlap_denominator = sqrt(irradiance_squared_integral_1300*irradiance_squared_integral_1800);
overlap_parameter = overlap_integral/overlap_denominator


omega_probe = 2*pi*3e8/(1300e-9);

n2_m2_per_W = 5e-18;
beta_2pa_cm_per_GW = 25;
beta_2pa_m_per_W = beta_2pa_cm_per_GW/1e11;
gamma_minv_Winv = omega_probe*overlap_parameter*n2_m2_per_W/(3e8*sqrt(mode_area_um2_1300*1e-12)*sqrt(mode_area_um2_1800*1e-12)) + 1i*overlap_parameter*beta_2pa_m_per_W/(sqrt(mode_area_um2_1300*1e-12)*sqrt(mode_area_um2_1800*1e-12))