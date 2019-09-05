function varargout = AudioTuningInspectionGUI_noise(varargin)
% AUDIOTUNINGINSPECTIONGUI_NOISE MATLAB code for AudioTuningInspectionGUI_noise.fig
%      AUDIOTUNINGINSPECTIONGUI_NOISE, by itself, creates a new AUDIOTUNINGINSPECTIONGUI_NOISE or raises the existing
%      singleton*.
%
%      H = AUDIOTUNINGINSPECTIONGUI_NOISE returns the handle to a new AUDIOTUNINGINSPECTIONGUI_NOISE or the handle to
%      the existing singleton*.
%
%      AUDIOTUNINGINSPECTIONGUI_NOISE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AUDIOTUNINGINSPECTIONGUI_NOISE.M with the given input arguments.
%
%      AUDIOTUNINGINSPECTIONGUI_NOISE('Property','Value',...) creates a new AUDIOTUNINGINSPECTIONGUI_NOISE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AudioTuningInspectionGUI_noise_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AudioTuningInspectionGUI_noise_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AudioTuningInspectionGUI_noise

% Last Modified by GUIDE v2.5 11-Jan-2019 11:27:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AudioTuningInspectionGUI_noise_OpeningFcn, ...
                   'gui_OutputFcn',  @AudioTuningInspectionGUI_noise_OutputFcn, ...
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

global data;


% --- Executes just before AudioTuningInspectionGUI_noise is made visible.
function AudioTuningInspectionGUI_noise_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AudioTuningInspectionGUI_noise (see VARARGIN)

% Choose default command line output for AudioTuningInspectionGUI_noise
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes AudioTuningInspectionGUI_noise wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = AudioTuningInspectionGUI_noise_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1
global data;

neuron = get(handles.listbox1,'Value');

cla(handles.axes1);
plot(handles.axes1,squeeze(data.eachTrial(neuron,:,:))','Color',[0.7 0.7 0.7]);
hold on;
plot(handles.axes1,data.responses(neuron,:),'Color',[0 0 0]);
xticks(handles.axes1,0:30:180);
xticklabels(handles.axes1,{'-1','0','1','2','3','4','5'});
xlabel(handles.axes1,'Time after stimulus onset (sec)');
ylabel(handles.axes1,'\Delta F/std(F)');


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global data;
% cla(handles.axes1);
% cla(handles.axes2);

[file folder] = uigetfile('/Users/Aaron/Documents');
data = load([folder '/' file]);

set(handles.text3,'String',file);

for i = 1:size(data.responses,1)
    list1{i} = num2str(i);
end
set(handles.listbox1,'String',list1);
