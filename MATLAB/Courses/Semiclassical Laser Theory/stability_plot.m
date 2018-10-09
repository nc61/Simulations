detuning_total = linspace(0,6,200);
Y1 = 1*ones(size(detuning_total));
detuning_y2 = detuning_total(detuning_total > 2);
Y2 = detuning_y2/2;
figure(1)
plot(detuning_total, Y1, detuning_y2, Y2), ylim([0,3]), xlabel('Detuning parameter \Delta'), ylabel('Field Intensity')