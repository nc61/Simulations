function [X, Y, beam_profile_flat, beam_profile_fit, fit_parameters] = fit_beam_profile_fcn(varargin)

if nargin == 2
    data_filename = file_selector('fileName', varargin{1});
    background_filename = file_selector('fileName', varargin{2});
    beam_profile_flat = read_camera_data_fcn(data_filename, background_filename);
else
    beam_profile_flat = read_camera_data_fcn();
end

x_range_scale = 4;
y_range_scale = 4;
    
pixel_pitch_um = 12.5;
res_x = 320;
res_y = 256;
 
[~, index_max] = max(beam_profile_flat);

[X,Y] = ndgrid(1:res_x, 1:res_y);

XY = zeros(length(beam_profile_flat),2);
XY(:,1) = X(:);
XY(:,2) = Y(:);

options = optimset('TolX', 1e-15);

fitfun = @(fit_parameters, XY)gaussian_2d(XY, fit_parameters(1), fit_parameters(2), fit_parameters(3), fit_parameters(4), fit_parameters(5), fit_parameters(6));

[peak_x_coordinate, peak_y_coordinate] = ind2sub([res_x, res_y],index_max);
fit_parameters = lsqcurvefit(fitfun, [1, 1, 1, 0, peak_x_coordinate, peak_y_coordinate], XY, beam_profile_flat, [0,0,0,0,0,0], [1.5, 320, 256, 360, 320,256], options);

trial_fit_max_value = fit_parameters(1);
trial_fit_x_waist = fit_parameters(2);
trial_fit_y_waist = fit_parameters(3);
trial_fit_angle = fit_parameters(4);
trial_fit_peak_x_coordinate = fit_parameters(5);
trial_fit_peak_y_coordinate = fit_parameters(6);

xrange = round(trial_fit_peak_x_coordinate + x_range_scale*trial_fit_x_waist*[-1,1]);
yrange = round(trial_fit_peak_y_coordinate + y_range_scale*trial_fit_y_waist*[-1,1]);

xrange = xrange(1):xrange(2);
yrange = yrange(1):yrange(2);
beam_profile_in_range = reshape(beam_profile_flat, 320, 256);
beam_profile_in_range = beam_profile_in_range(xrange, yrange);
beam_profile_in_range_flat = beam_profile_in_range(:);

[X,Y] = ndgrid(xrange, yrange);

XY = zeros(length(beam_profile_in_range_flat),2);
XY(:,1) = X(:);
XY(:,2) = Y(:);

fit_parameters = lsqcurvefit(fitfun, [trial_fit_max_value, trial_fit_x_waist, trial_fit_y_waist, trial_fit_angle, trial_fit_peak_x_coordinate, trial_fit_peak_y_coordinate], XY, beam_profile_in_range_flat, [0,0,0,0,0,0], [1.5, 320, 256, 360, 320,256], options);

beam_profile_fit = reshape(fitfun(fit_parameters, XY), [length(xrange), length(yrange)]);
beam_profile_flat = beam_profile_in_range;
X = X*pixel_pitch_um;
Y = Y*pixel_pitch_um;


end

