filenames = {'p0.raw', 'p2.raw', 'p4.raw', 'p6.raw', 'p8.raw'};
backgrounds = {'p1.raw', 'p3.raw', 'p5.raw', 'p7.raw', 'p9.raw'};

x_waist = zeros(size(filenames));
y_waist = zeros(size(filenames));
X = cell(1,length(filenames));
Y = cell(1,length(filenames));
beam_profile = cell(1,length(filenames));
beam_profile_fit = cell(1,length(filenames));
fit_parameters = cell(1,length(filenames));

for ind = 1:length(filenames)
    
    [X{ind}, Y{ind}, beam_profile{ind}, beam_profile_fit{ind}, fit_parameters{ind}] = fit_beam_profile_fcn(filenames{ind}, backgrounds{ind});
    fit_params = fit_parameters{ind};
    x_waist(ind) = fit_params(2)*12.5e-3;
    y_waist(ind) = fit_params(3)*12.5e-3;
    

end

index_to_plot = 2;

figure(1)
surf(X{index_to_plot},Y{index_to_plot}, beam_profile{index_to_plot})
fit_params = fit_parameters{index_to_plot};
colormap jet, shading interp, caxis([0 fit_params(1)])

figure(2)
surf(X{index_to_plot},Y{index_to_plot}, beam_profile_fit{index_to_plot})
colormap jet, shading interp
% 
% figure(3)
% plot([10], x_waist, 'ro');
% hold on
% plot([10], y_waist, 'bo');
% hold off
