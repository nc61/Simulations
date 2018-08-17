function full_filepath = file_selector(varargin)

p = inputParser;

addParameter(p, 'fileName', '');
addParameter(p, 'selectDirectory', 0);
addParameter(p, 'multiSelect', 'off');
parse(p, varargin{:})

filename = p.Results.fileName;
select_multiple = p.Results.multiSelect;
select_directory = p.Results.selectDirectory;

call_stack = dbstack;
calling_function = call_stack(end).name;
[file_selector_directory, ~, ~] = fileparts(mfilename('fullpath'));
file_selector_userdata_path = fullfile(file_selector_directory, 'file_selector_userdata');

if exist(file_selector_userdata_path, 'file') == 2
    userdata_text = fileread(file_selector_userdata_path);
else
    userdata_text = '';
end

if ~isempty(userdata_text)
    userdata_lines = regexp(userdata_text, '\n', 'split');
    userdata_lines_split = regexp(userdata_lines, '\t', 'split');
    function_names = cellfun(@(line) line(1), userdata_lines_split);
    function_line = find(strcmp(function_names, calling_function));
else
    userdata_lines = {};
    userdata_lines_split = {};
    function_line = '';
end

if ~isempty(function_line)
    function_saved_path = fileparts(userdata_lines_split{function_line}{2});
else
    function_saved_path = '';
    function_line = length(userdata_lines) + 1;
end

if ~isempty(filename)
    if ~select_directory
        full_filepath = fullfile(function_saved_path, filename);
        if exist(full_filepath, 'file') ~= 2
            [~, pathname] = uigetfile('*.*');
            full_filepath = fullfile(pathname, filename);
            if exist(full_filepath, 'file') ~= 2
                error('The file %s does not exist in the directory %s', filename, pathname);
            end
        else
            pathname = function_saved_path;
        end
    end
    
else
    if strcmpi(select_multiple, 'on')
         [filename, pathname] = uigetfile({sprintf('%s\\*.*', function_saved_path)}, 'Multiselect', 'on');
    else
        [filename, pathname] = uigetfile({sprintf('%s\\*.*', function_saved_path)});
    end
        full_filepath = fullfile(pathname, filename);
end

userdata_lines_split{function_line}{1} = calling_function;
userdata_lines_split{function_line}{2} = pathname;

fid_write = fopen(file_selector_userdata_path, 'w');

for ind = 1:length(userdata_lines_split)-1
    line_to_write = sprintf('%s\t%s\n', userdata_lines_split{ind}{1}, userdata_lines_split{ind}{2});
    fwrite(fid_write, line_to_write);
end

line_to_write = sprintf('%s\t%s', userdata_lines_split{length(userdata_lines_split)}{1}, userdata_lines_split{length(userdata_lines_split)}{2});
fwrite(fid_write, line_to_write);
fclose(fid_write);


end


