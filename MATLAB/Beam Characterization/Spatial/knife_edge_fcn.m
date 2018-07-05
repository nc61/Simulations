function [waist_HW_1_over_e_max_mm, x_mm, signal_arb] = knife_edge_fcn(varargin)

if nargin == 1
    filename = uigetfile();
    plot_data_flag = varargin{1};
else
    filename = varargin{1};
    plot_data_flag = varargin{2};
end

scan_data = tdfread(filename);
stage_position = scan_data.Position;
signal_arb = scan_data.x1;

if (signal_arb(1) < signal_arb(2))
    fitfun = @(p, xdata)p(1) + p(2)/2*(1 - erf((xdata - p(3))/p(4)));
else
    fitfun = @(p, xdata)p(1) + p(2)/2*(1 + erf((xdata - p(3))/p(4)));
end

maximum = 20;
midpoint_absorption_value = (signal_arb(1) + signal_arb(end))/2;
[~, stage_midpoint_index] = min(abs(signal_arb - midpoint_absorption_value));
center = stage_position(stage_midpoint_index);
waist_HW_1_over_e_max_mm = (stage_position(end) - stage_position(1))/5;


fit_parameters = lsqcurvefit(fitfun, [offset, maximum, center, waist_HW_1_over_e_max_mm], stage_position, signal_arb);

xfit = linspace(stage_position(1), stage_position(end), 1000);

plot(stage_position, signal_arb, 'bo')
hold on
plot(xfit, fitfun(fit_parameters,xfit), 'r')
hold off

waist_HW_1_over_e_max_mm = fit_parameters(4)