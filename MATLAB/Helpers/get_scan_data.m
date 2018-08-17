function [stage_position, signal] = get_scan_data(filename, varargin)

p = inputParser;
addOptional(p, 'channel', 1);
parse(p, varargin{:});

channel = p.Results.channel;

scan_data = dlmread(filename, '\t', 1,0);
stage_position = scan_data(:,1);
signal = scan_data(:, 1 + channel);


