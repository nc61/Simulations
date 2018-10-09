function [waist_HW_1_over_e_max_mm, x_mm, signal] = knife_edge_fcn(varargin)

p = inputParser;
addOptional(p, 'fileName', '');
addParameter(p, 'inputXData', []);
addParameter(p, 'inputYData', []);
addParameter(p, 'average', 'off')
addParameter(p, 'figureNum', 1)

parse(p, varargin{:});
filename = p.Results.fileName;
stage_position = p.Results.inputXData;
signal = p.Results.inputYData;
average_flag = p.Results.average;
figure_num = p.Results.figureNum;


if isempty(signal)
    
    if strcmpi(average_flag, 'on')
        [stage_position, signal] = average_scans_fcn();
    else
        full_filepath = file_selector('fileName', filename);
        [stage_position, signal] = get_scan_data(full_filepath);
    end    

end

if (signal(1) < signal(2))
    fitfun = @(p, xdata)p(1) + p(2)/2*(1 - erf(sqrt(2)*(xdata - p(3))/p(4)));
else
    fitfun = @(p, xdata)p(1) + p(2)/2*(1 + erf(sqrt(2)*(xdata - p(3))/p(4)));
end

maximum = max(signal);
midpoint_absorption_value = (signal(1) + signal(end))/2;
[~, stage_midpoint_index] = min(abs(signal - midpoint_absorption_value));
center = stage_position(stage_midpoint_index);
waist_HW_1_over_e_max_mm = (stage_position(end) - stage_position(1))/5;
offset = min(abs(signal));

fit_parameters = lsqcurvefit(fitfun, [offset, maximum, center, waist_HW_1_over_e_max_mm], stage_position, signal);

xfit = linspace(stage_position(1), stage_position(end), 1000);

figure(figure_num)
plot(stage_position, signal, 'bo')
hold on
plot(xfit, fitfun(fit_parameters,xfit), 'r')
hold off
xlabel('knife position (mm)'), ylabel('Energy (arb)'), title('Knife edge scan')

waist_HW_1_over_e_max_mm = fit_parameters(4);
x_mm = stage_position;