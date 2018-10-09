function [mode_area_1_um2, mode_area_2_um2, overlap_parameter] = mode_area_and_overlap(filename_1, mode_number_1, polarization_1, filename_2, mode_number_2, polarization_2)

mode_1 = load(filename_1, sprintf('mode%d_E%s', mode_number_1, polarization_1));
mode_1_fields = fieldnames(mode_1);
mode_1 = getfield(mode_1, mode_1_fields{1});
load(filename_1, 'x', 'y', 'z');

mode_2 = load(filename_2, sprintf('mode%d_E%s', mode_number_2, polarization_2));
mode_2_fields = fieldnames(mode_2);
mode_2 = getfield(mode_2, mode_2_fields{1});

mode_irradiance_1 = abs(mode_1).^2;
mode_irradiance_squared_1 = abs(mode_1).^4;
irradiance_integral_1 = trapz(z, trapz(y, mode_irradiance_1, 2), 3);
irradiance_squared_integral_1 = trapz(z, trapz(y, mode_irradiance_squared_1, 2), 3);
mode_area_1_m2 = irradiance_integral_1^2/irradiance_squared_integral_1;
mode_area_1_um2 = mode_area_1_m2*1e12;

mode_irradiance_2 = abs(mode_2).^2;
mode_irradiance_squared_2 = abs(mode_2).^4;
irradiance_integral_2 = trapz(z, trapz(y, mode_irradiance_2, 2), 3);
irradiance_squared_integral_2 = trapz(z, trapz(y, mode_irradiance_squared_2, 2), 3);
mode_area_2_m2 = irradiance_integral_2^2/irradiance_squared_integral_2;
mode_area_2_um2 = mode_area_2_m2*1e12;

overlap_integral = trapz(z, trapz(y, mode_irradiance_1.*mode_irradiance_2, 2), 3);
overlap_denominator = sqrt(irradiance_squared_integral_1*irradiance_squared_integral_2);
overlap_parameter = overlap_integral/overlap_denominator;


end

