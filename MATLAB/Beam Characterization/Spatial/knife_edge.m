scan_data = tdfread('H:\Home\NLO\Shared\Group Members\Nick\Measurements 071918\scan7');
stage_position = scan_data.Position;
absorp = scan_data.x1;

if (absorp(1) < absorp(2))
    fitfun = @(p, xdata)p(1) + p(2)/2*(1 - erf((xdata - p(3))/p(4)));
else
    fitfun = @(p, xdata)p(1) + p(2)/2*(1 + erf((xdata - p(3))/p(4)));
end

maximum = 5;
midpoint_absorption_value = (absorp(1) + absorp(end))/2;
[~, stage_midpoint_index] = min(abs(absorp - midpoint_absorption_value));
center = stage_position(stage_midpoint_index);
waist_HW_1_over_e_max = (stage_position(end) - stage_position(1))/5;
offset = 0;


fit_parameters = lsqcurvefit(fitfun, [offset, maximum, center, waist_HW_1_over_e_max], stage_position, absorp);

xfit = linspace(stage_position(1), stage_position(end), 1000);

plot(stage_position, absorp, 'bo')
hold on
plot(xfit, fitfun(fit_parameters,xfit), 'r')
hold off

waist_HW_1_over_e_max = fit_parameters(4)