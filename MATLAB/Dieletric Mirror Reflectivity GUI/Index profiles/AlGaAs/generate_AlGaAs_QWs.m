function generate_AlGaAs_QWs(min_wavelength_um, max_wavelength_um, well_comp, well_width_nm, barrier_comp, barrier_width_nm)

lambda_m = min_wavelength_um*1e-6:1e-9:max_wavelength_um*1e-6;

n_well = index2(lambda_m, well_comp);
n_barrier = index2(lambda_m, barrier_comp);

n_average = sqrt((well_width_nm*n_well.^2 + barrier_width_nm*n_barrier.^2)./(well_width_nm + barrier_width_nm));
fname = "QW_"+num2str(well_width_nm)+ "nm"+"(x="+num2str(well_comp)+")_barrier_"+num2str(barrier_width_nm)+"nm(x="+num2str(barrier_comp)+").txt";

fid = fopen(fname, 'w+');
k = 0;
format long
for jnd = 1:length(lambda_m)
    line_to_write = sprintf('%s\t%s\t%s\r\n', num2str(lambda_m(jnd)*1e6), n_average(jnd), num2str(k));
    fwrite(fid, line_to_write);
end
fclose(fid);

