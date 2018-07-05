fname = 'found zero delay';
probe = 3.42*4;

scan_data = tdfread(fname);
    stage_position = scan_data.Position;
    absorp = scan_data.x1;
    
    plot(stage_position, absorp);
    
    [~, zero_delay_index] = max(absorp);
    delay = 2*(stage_position - stage_position(zero_delay_index))/10*2.54/100/2.998e8*1e15; % [fs]
    absorption = (absorp - mean(absorp(1:20)))/probe_array(ind);
    

figure(2)
plot(delay, absorption, 'r*')
xlabel('delay (fs)'), ylabel('\DeltaE/E'), xlim([min(delay) max(delay)]), title('Pump = 1200nm, Probe = 2600nm')