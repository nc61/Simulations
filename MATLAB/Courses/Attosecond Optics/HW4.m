A = 1.4580;
B = 0.00354;

n = @(lambda_um) A + B./lambda_um.^2;
c = 3e8;

n_800 = n(0.8);
vp_800 = c/n_800

n_1600 = n(1.6);
vp_1600 = c/n_1600

vg_800 = c/(n(0.8) + 2*B/0.8^2)
vg_1600 = c/(n(1.6) + 2*B/1.6^2)

b2_800 = 0.8e-6^3/(2*pi*c^2)*6*B*(1e-6)^2/(0.8e-6^4)*1e15^2/1000
b2_1600 = 1.6e-6^3/(2*pi*c^2)*6*B*(1e-6)^2/(1.6e-6^4)*1e15^2/1000

lambda = 0.4:0.01:2
n_l = n(lambda);
figure(1)
plot(lambda, n_l)

d_lambda = lambda(2) - lambda(1);
lambda_deriv = diff(n_l)/d_lambda;
figure(2)
plot(lambda_deriv)

lambda_2nd_deriv = diff(lambda_deriv)/d_lambda;

ng = n_l(lambda == 0.8) - 0.8*lambda_deriv(lambda == 0.8)
c/ng

beta = 1e-6*0.8^3/(2*pi*3e8^2)*lambda_2nd_deriv(lambda == 0.8)*1e15^2/1e3