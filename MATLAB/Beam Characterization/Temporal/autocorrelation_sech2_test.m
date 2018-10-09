fitfun = @(p, xdata) p(3) + p(4)*4*csch((xdata - p(2))/p(1)).^2.*(-1 + ((xdata - p(2))/p(1)).*coth((xdata - p(2))/p(1)));

x = linspace(-10,10,200);
fit_params = [1 0 0.2 1];

func_results = fitfun(fit_params, x);
plot(x,func_results)