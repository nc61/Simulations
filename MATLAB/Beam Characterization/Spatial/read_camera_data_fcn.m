function data = read_camera_data_fcn(varargin)

if nargin == 1
    filename = varargin{1};
    background_filename = '';
elseif nargin == 2
    filename = varargin{1};
    background_filename = varargin{2};
else
    error('No file to read');
end

fid = fopen(filename);
data = fread(fid, 'uint16=>double');
fclose(fid);


if ~isempty(background_filename)
    fid_background = fopen(background_filename);
    data_background = fread(fid_background, 'uint16=>double');
    fclose(fid_background);
    
    data = data - data_background;
    
end

data = data/max(max(data));


end

