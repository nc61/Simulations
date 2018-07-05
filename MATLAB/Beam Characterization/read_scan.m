
fname_array = ["scan1", "scan2", "scan3", "scan4", "scan5", "scan6", "scan7", "scan8"];
probe_array = [7, 3, 3.2, 7.2, 4.3, 4.8, 2.3, 3.2];

delay = cell(length(fname_array),1);
absorption = cell(length(fname_array),1);

for ind = 1:length(fname_array)
    scan_data = tdfread(fname_array(ind));
    stage_position = scan_data.Position;
    absorp = scan_data.x2;
    
    plot(stage_position, absorp);
    
    [~, zero_delay_index] = min(absorp);
    delay{ind} = -2*(stage_position - stage_position(zero_delay_index))/1000/2.998e8*1e15; % [fs]
    absorption{ind} = (absorp - mean(absorp(1:20)))/probe_array(ind);
    
end

si_scan_indices = [1,4,6,8];
wavelengths = [1300, 1340, 1400, 1440];

signal_strength = zeros(size(si_scan_indices));
for jnd = 1:length(si_scan_indices)
    signal_strength(jnd) = abs(min(absorption{si_scan_indices(jnd)}));
end

figure(1)
plot(1240./wavelengths, signal_strength, 'o')


figure(2)
plot(delay{8}, absorption{8})