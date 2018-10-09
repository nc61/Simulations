function signal = sech2_fit_function(p, xdata)

signal = p(3) + p(4)*4*csch((xdata - p(2))/p(1)).^2.*(-1 + ((xdata - p(2))/p(1)).*coth((xdata - p(2))/p(1)));
signal(xdata == p(2)) = p(4)*4/3;

end

