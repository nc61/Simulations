
c = 3e8;
n = transpose(real(neff));

f_interp = linspace(min(f), max(f), 1000);
n_interp = interp1(f, n, f_interp);

dn_df = [0 diff(n_interp)]/(f_interp(2) - f_interp(1));
vg_f = c./(n_interp + f_interp.*dn_df);
beta_2_f = 1/(2*pi)*[0 diff(1./vg_f)]/(f_interp(2) - f_interp(1))*(1e15)^2/1000;

plot(f_interp, beta_2_f), ylim([0 2000]);