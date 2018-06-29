L = 50;
N = 2^10;
x = (L/N)*(0:N-1);
kx = 2*pi/L*[(0:N/2-1) (-N/2:-1)];

Z = 20;
steps = 1000;
dZ = Z/steps;
nl = 1;
alpha = 0;


phi = exp(-0.04*(x - 0.5*L).^2)
V = 0;

D = -0.5*(1i*kx.^2 + alpha)*0.5*dZ;

phi_out = zeros(steps, N);

for ind = 1:steps
    lin_step_1 = ifft(fft(phi).*exp(D));
    %nonlinear_step = lin_step_1.*exp((1i*V + 1i.*nl*abs(lin_step_1).^2)*dZ);
    nonlinear_step = lin_step_1;
    lin_step_2 = ifft(fft(nonlinear_step).*exp(D));
    phi = lin_step_2;
    phi_out(ind, :) = phi;
end

P = zeros(1,steps);
%Power
dX = L/N;
for i = 1:steps
    P(i) = dX*sum(abs(phi_out(i,:)).^2);
end
 
%Hamiltonian
%Numerical differentiation
dphi_dx = zeros(steps,N);
for i = 1:steps
    dphi_dx(i,1) = ((phi_out(i,2))-phi_out(i,1))/dX;
    dphi_dx(i,2) = (phi_out(i,3)-phi_out(i,1))/(2*dX);
end
for i = 1:steps
    dphi_dx(i,end) = (phi_out(i,end)-(phi_out(i,end-1)))/dX;
    dphi_dx(i,end-1) = (phi_out(i,end)-(phi_out(i,end-2)))/(2*dX);
end

for j = 3:N-2
for i = 1:steps
    dphi_dx(i,j) = (-phi_out(i,j+2)+8*phi_out(i,j+1)-8*phi_out(i,j-1)+phi_out(i,j-2))/(12*dX);
end
end
 
H = zeros(1, steps);
%Calculation of Hamiltonian
for i = 1:steps
    H(i) = dX*sum(abs(dphi_dx(i,:)).^2-nl*abs(phi_out(i,:)).^4-2*V.*abs(phi_out(i,:)).^2);
end


z = linspace(0,Z,steps);
figure(1)
surf(x,z,abs(phi_out).^2)
shading interp
view(2)
colormap('jet')
colorbar
xlabel('T (a.u.)')
ylabel('Distance (a.u.)')
xlim([min(x) max(x)])

%Power conservation
figure(2)
plot(z,P,'Linewidth',2)
xlabel('Distance (a.u.)')
ylabel('P')
ylim([P(1)-1 P(1)+1])
 
%Hamiltonian
figure(3)
plot(z,H,'Linewidth',2)
xlabel('Distance (a.u.)')
ylabel('H')
ylim([H(1)-10 H(1)+10])

 
