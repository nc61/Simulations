v0 = 1/4;
x = linspace(-8,8,200);
y = linspace(-8,8,200);
[X,Y] = ndgrid(x,y);
f1 = 2*cos(2*pi*v0*X) + 2*cos(2*pi*v0*Y) + 4*cos(2*pi*v0/sqrt(2)*X).*cos(2*pi*v0/sqrt(2)*Y);

figure(1)
surf(X,Y,abs(f1).^2), shading interp, xlabel('x/\lambda'), ylabel('y/\lambda')

z = 100;
f2 = 2*exp(1i*pi*z*v0^2)*(cos(2*pi*v0*X) + cos(2*pi*v0*Y)) + 4*exp(1i*pi*z*v0^2)*(cos(2*pi*v0/sqrt(2)*X).*cos(2*pi*v0/sqrt(2)*Y));

figure(2)
surf(X,Y,abs(f2).^2), shading interp, xlabel('x/\lambda'), ylabel('y/\lambda')
