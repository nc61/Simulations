function [position] = read_tab_delimited_file(filename, headerlines)

text_data = fileread(filename);
lines = regexp(text_data, '\r\n', 'split');

if isempty(lines{end})
    lines = lines(1:end-1);
end

lines_without_header = lines(1 + headerlines: end); 
split_lines = regexp(lines_without_header, '\t', 'split');


position = str2double(cellfun(@(element) element(1), split_lines));


end

