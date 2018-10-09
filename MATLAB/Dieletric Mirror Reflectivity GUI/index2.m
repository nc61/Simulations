function n = index2(lambda, x)

h = 6.626e-34;
c = 2.998e8;
q = 1.6e-19;

E = @(lambda) h*c./(lambda*q);
Ep = E(lambda);


A = @(x) 5.6684 + 10.464.*x + 1.450*exp(3.5584.*x);
B = @(x) 9.1813 - 4.5059.*x - 1.4304.*x.^2 + 2.2338.*x.^3;
e = @(x)2.1582 + 0.80331.*x - 0.0911.*x.^2 - 2.6906.*x.^3;


E0 = @(x) 1.422 + 1.2475.*x;

delta0 = @(x) 0.34 - 0.02.*x;

E1 = @(x) 2.926 + 0.6717.*x - 0.3242.*x.^2 + 0.6172.*x.^3;


f = @(y)(2 - (1 + y).^(0.5) - (1 - y).^(0.5))./(y.^2);
chi = @(Et) Ep/Et;
g = (Ep./E1(x)).^2.*log((1 - (Ep.^2./E1(x).^2))/(1 - Ep.^2./6.4^2));

eps_0 = e(x);
eps_1 = A(x).*(E0(x)).^(-1.5).*(f(chi(E0(x))) + (1/2)*(E0(x)./(E0(x) + delta0(x))).^(1.5).*f(chi(E0(x) + delta0(x))));
eps_2 = -B(x).*(chi(E1(x))).^(-2).*log((1 - chi(E1(x)).^2)./(1 - chi(6.4).^2));

eps = eps_0 + eps_1 + eps_2;
n = real(sqrt(eps));
end



