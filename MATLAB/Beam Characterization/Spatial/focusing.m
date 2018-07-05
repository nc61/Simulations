positions = [0 5 10 15 20 25];
waists = zeros(size(positions));

filenames = {'e_0', 'e_5', 'e_10', 'e_15', 'e_20', 'e_25'};

for ind = 1:length(filenames)

scan_data = tdfread(filenames{ind});
stage_position = scan_data.Position;
absorp = scan_data.x1;

if (absorp(1) < absorp(2))
    fitfun = @(p, xdata)p(1) + p(2)/2*(1 - erf((xdata - p(3))/p(4)));
else
    fitfun = @(p, xdata)p(1) + p(2)/2*(1 + erf((xdata - p(3))/p(4)));
end

maximum = 8;
midpoint_absorption_value = (absorp(1) + absorp(end))/2;
[~, stage_midpoint_index] = min(abs(absorp - midpoint_absorption_value));
center = stage_position(stage_midpoint_index);
center = 4.5
waist = (stage_position(end) - stage_position(1))/5;
waist = 0.6
offset = 0;


fit_parameters = lsqcurvefit(fitfun, [offset, maximum, center, waist], stage_position, absorp);

xfit = linspace(stage_position(1), stage_position(end), 1000);

figure(ind)
plot(stage_position, absorp, 'bo')
xlabel('position (mm)'), ylabel('signal (a.u.)'), title(sprintf('knife edge at x = %dmm', positions(ind)))
hold on
plot(xfit, fitfun(fit_parameters,xfit), 'r')
hold off

waist = fit_parameters(4);

waists(ind) = waist;

end

figure(length(filenames) + 1);
plot(positions, waists, '*');