%%%%% Problem 5

c = 3e8;
eps0 = 8.85e-12;
lambda = 1.53*1e-6; %m
k0 = 2*pi/lambda;
D = 5*1e-12*1e-3*1e9; %s/(m m)
tau_0 = 1e-12;
Tr = 4e-15;

beta_0 = -(lambda)^2/(2*pi*c)*D
dw = 2*pi*c*(1/(lambda + 10e-9) - 1/(lambda))

z = dw*(-15/8)*tau_0^4/(Tr)*1/abs(beta_0)

d2 = 2000;
alpha = -log(10^(-2/10))/1000;
zeff = (1 - exp(-alpha*d2))/alpha
dlambda = -lambda^2/(2*pi*c)*8/15*(Tr/tau_0^4)*abs(beta_0)*zeff