c = 3e8;

lambda = linspace(1e-6, 2e-6, 5000);

f = linspace(c/2e-6, c/1e-6, 1000); 

refractive_index = index2(lambda, 0.16);
refractive_index_f = index2(c./f, 0.16);


omega = 2*pi*3e8./lambda;
k0 = 2*pi./lambda;
dn_dlambda = [0 diff(refractive_index)]/(lambda(2) - lambda(1));
dn_df = [0 diff(refractive_index_f)]/(f(2) - f(1));
c = 3e8;

vg = c./(refractive_index - lambda.*dn_dlambda);
vg_f = c./(refractive_index_f + f.*dn_df);
beta_2_f = 1/(2*pi)*[0 diff(1./vg_f)]/(f(2) - f(1))*(1e15)^2/1000;

beta_2 = -lambda.^2/(2*pi*c).*[0 diff(1./vg)]/(lambda(2) - lambda(1))*(1e15)^2/1000;
figure(1)
plot(lambda, vg)
hold on
plot(c./f, vg_f, 'r')
hold off

figure(2)
plot(lambda, beta_2), xlim([1.01e-6 2e-6]), ylim([0 5000])
hold on
plot(c./f, beta_2_f, 'r')
hold off