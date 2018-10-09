syms w0 w gamma z

f1 = 1i*pi*(z - w0)/(((z - w0)^2 + gamma^2));
f1 = subs(f1, z, w)

f2 = 1i*2*pi*(z - w0)/((z - (w0 - 1i*gamma))*(z - w));
f2 = subs(f2, z, w0 + 1i*gamma)