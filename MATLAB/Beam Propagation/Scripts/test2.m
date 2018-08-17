filenames = {'without_filtering.raw'};
backgrounds = {'without_filtering1.raw'};

x_waist = zeros(size(filenames));
y_waist = zeros(size(filenames));

for ind = 1:length(filenames)
   
    [X{ind}, Y{ind}, beam_profile{ind}, beam_profile_fit{ind}, fit_parameters] = fit_beam_profile_fcn(filenames{ind}, backgrounds{ind});
    x_waist(ind) = fit_parameters(2)*12.5e-3;
    y_waist(ind) = fit_parameters(3)*12.5e-3;
   

end

    figure(1)
    surf(X,Y, beam_profile)
    colormap jet, shading interp
    
    figure(2)
    surf(X,Y, beam_profile_fit)
    colormap jet, shading interp
    