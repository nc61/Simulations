function F_out = F(n, sign, me_z, mh_z, d, Ep_1, Ep_2)
hbar = 6.626e-34/(2*pi);


if sign=='+'
    F_out = ((4*(n+1)*n)./(1+2*n)).^2.*...
         ((Ep_1./((Ep_1 - (1 + 2*n)*pi^2*hbar^2/(2*me_z*d^2)).*(Ep_1 + (1 + 2*n)*pi^2*hbar^2/(2*mh_z*d^2))))+...
         (Ep_2./((Ep_2 - (1 + 2*n)*pi^2*hbar^2/(2*me_z*d^2)).*(Ep_2 + (1 + 2*n)*pi^2*hbar^2/(2*mh_z*d^2))))).^2;
elseif sign=='-'
    F_out = ((4*(n-1)*n)/(1-2*n)).^2.*...
        ((Ep_1./((Ep_1 - (1 - 2*n)*pi^2*hbar^2/(2*me_z*d^2)).*(Ep_1 + (1 - 2*n)*pi^2*hbar^2/(2*mh_z*d^2))))+...
        (Ep_2./((Ep_2 - (1 - 2*n)*pi^2*hbar^2/(2*me_z*d^2)).*(Ep_2 + (1 - 2*n)*pi^2*hbar^2/(2*mh_z*d^2))))).^2;
else
    F_out=0;
end

