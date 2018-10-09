function varargout = reflectivity_gui(varargin)
% REFLECTIVITY_GUI MATLAB code for reflectivity_gui.fig
%      REFLECTIVITY_GUI, by itself, creates a new REFLECTIVITY_GUI or raises the existing
%      singleton*.
%
%      H = REFLECTIVITY_GUI returns the handle to a new REFLECTIVITY_GUI or the handle to
%      the existing singleton*.
%
%      REFLECTIVITY_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in REFLECTIVITY_GUI.M with the given input arguments.
%
%      REFLECTIVITY_GUI('Property','Value',...) creates a new REFLECTIVITY_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before reflectivity_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to reflectivity_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help reflectivity_gui

% Last Modified by GUIDE v2.5 21-Aug-2015 13:18:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @reflectivity_gui_OpeningFcn, ...
    'gui_OutputFcn',  @reflectivity_gui_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before reflectivity_gui is made visible.
function reflectivity_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to reflectivity_gui (see VARARGIN)

% Choose default command line output for reflectivity_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes reflectivity_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);

set(handles.wv_min_edit, 'String', sprintf('%d', 850));
set(handles.wv_max_edit, 'String', sprintf('%d', 1150));
set(handles.Al_scale_edit, 'String', sprintf('%d', 100));
set(handles.Ga_scale_edit, 'String', sprintf('%d', 100));

handles.measured_data = [];
handles.layers = [];
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = reflectivity_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in lfile_pushbutton.
function lfile_pushbutton_Callback(hObject, eventdata, handles)
[fname, pname] = uigetfile('*.*', 'Choose a text file...');
full_filename = fullfile(pname, fname);
set(handles.lfile_edit, 'String', sprintf('%s%s', pname, fname));

wv_min = str2double(get(handles.wv_min_edit, 'String'))*1e-9;
wv_max = str2double(get(handles.wv_max_edit, 'String'))*1e-9;

[wv, R, layers] = reflectivity_fcn([wv_min wv_max], 101, 1, full_filename);
axes(handles.reflectivity_axes);
plot(wv, R);
xlim([wv_min wv_max]);

handles.layer_data = layers;
guidata(hObject, handles);

function lfile_edit_Callback(hObject, eventdata, handles)



function Al_scale_edit_Callback(hObject, eventdata, handles)

Al_scale = str2double(get(handles.Al_scale_edit, 'String'));

if (Al_scale < 0)
    set(handles.Al_scale_edit, 'String', sprintf('%d', 0));
end

update_plot(hObject, handles);

function Ga_scale_edit_Callback(hObject, eventdata, handles)
Ga_scale = str2double(get(handles.Ga_scale_edit, 'String'));

if (Ga_scale < 0)
    set(handles.Ga_scale_edit, 'String', sprintf('%d', 0));
end

update_plot(hObject, handles);



function update_plot(hObject, handles)

layers = handles.layer_data;

if (~isempty(layers))
    Ga_scale = str2double(get(handles.Ga_scale_edit, 'String'));
    Al_scale = str2double(get(handles.Al_scale_edit, 'String'));
    
    layers(:, 5) = (layers(:,3).*(Ga_scale/100) + layers(:,4).*(Al_scale/100)).*layers(:, 5);
    
    wv_min = str2double(get(handles.wv_min_edit, 'String'))*1e-9;
    wv_max = str2double(get(handles.wv_max_edit, 'String'))*1e-9;
    
    [wv, R, ~] = reflectivity_fcn([wv_min wv_max], 101, 0, layers);
    axes(handles.reflectivity_axes);
    fclose all;
    
    plot(wv, R);
end

if (~isempty(handles.measured_data))
    measured_data_shaped = handles.measured_data;
    hold on
    plot(1e-9*measured_data_shaped(:,1), measured_data_shaped(:,2)/100, 'r');
    legend('calculated', 'measured');
    xlabel('wavelength (nm)'), ylabel('R');
    hold off
    handles.measured_data = measured_data_shaped;
    guidata(hObject, handles);
end
xlim([wv_min wv_max]);


function wv_min_edit_Callback(hObject, eventdata, handles)
update_plot(hObject, handles);



function wv_max_edit_Callback(hObject, eventdata, handles)
update_plot(hObject, handles)



function lfile_ref_edit_Callback(hObject, eventdata, handles)



% --- Executes on button press in lfile_ref_pushbutton.
function lfile_ref_pushbutton_Callback(hObject, eventdata, handles)
[fname, pname] = uigetfile('*.*', 'Choose a text file...');
full_filename = fullfile(pname, fname);
set(handles.lfile_ref_edit, 'String', sprintf('%s%s', pname, fname));

fid = fopen(full_filename);
measured_data = fscanf(fid, '%f');
measured_data_shaped = reshape(measured_data, 2, length(measured_data)/2)';
fclose(fid);

hold on
plot(1e-9*measured_data_shaped(:,1), measured_data_shaped(:,2)/100, 'r');
hold off

handles.measured_data = measured_data_shaped;
guidata(hObject, handles);


% --- Executes on slider movement.
function Al_scale_slider_Callback(hObject, eventdata, handles)
scale = get(handles.Al_scale_slider, 'Value');
set(handles.Al_scale_edit, 'String', sprintf('%.2f', scale));
update_plot(hObject, handles);


% --- Executes on slider movement.
function Ga_scale_slider_Callback(hObject, eventdata, handles)
scale = get(handles.Ga_scale_slider, 'Value');
set(handles.Ga_scale_edit, 'String', sprintf('%.2f', scale));
update_plot(hObject, handles);
