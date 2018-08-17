filename = file_selector();
[stage_position, signal] = get_scan_data(filename);

plot(2.54*stage_position, signal, 'r')
xlabel('stage_position (mm)'), ylabel('signal (arb)')

