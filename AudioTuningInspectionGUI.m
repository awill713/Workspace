function varargout = AudioTuningInspectionGUI_noAttn(varargin)
% AUDIOTUNINGINSPECTIONGUI MATLAB code for AudioTuningInspectionGUI.fig
%      AUDIOTUNINGINSPECTIONGUI, by itself, creates a new AUDIOTUNINGINSPECTIONGUI or raises the existing
%      singleton*.
%
%      H = AUDIOTUNINGINSPECTIONGUI returns the handle to a new AUDIOTUNINGINSPECTIONGUI or the handle to
%      the existing singleton*.
%
%      AUDIOTUNINGINSPECTIONGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AUDIOTUNINGINSPECTIONGUI.M with the given input arguments.
%
%      AUDIOTUNINGINSPECTIONGUI('Property','Value',...) creates a new AUDIOTUNINGINSPECTIONGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AudioTuningInspectionGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AudioTuningInspectionGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AudioTuningInspectionGUI

% Last Modified by GUIDE v2.5 28-May-2018 14:33:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AudioTuningInspectionGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @AudioTuningInspectionGUI_OutputFcn, ...
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




% --- Executes just before AudioTuningInspectionGUI is made visible.
function AudioTuningInspectionGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AudioTuningInspectionGUI (see VARARGIN)

% Choose default command line output for AudioTuningInspectionGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes AudioTuningInspectionGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = AudioTuningInspectionGUI_OutputFcn(hObject, eventdata, handles) 
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
attn = get(handles.listbox3,'Value');
tone = get(handles.listbox2,'Value');

cla(handles.axes1);
for a = 1:length(data.attenuations)
    tCurve = squeeze(data.tuningCurves(neuron,a,:));

    semilogx(handles.axes1,data.frequencies,tCurve);
    xlabel(handles.axes1,'Frequency (Hz)');
    ylabel(handles.axes1,'\Delta F/std(F)');
    hold(handles.axes1,'on');
end
legend(handles.axes1,data.attnString);

% ylim(handles.axes1,[-0.2 1]);
% set(handles.axes1,'xlim',[3 5],'xticks',[3 4 5],'xticklabels',{1 2 3},'xlabel','Frequency (kHz)');
% xticks(handles.axes1,[3 4 5]
% xticklabels(handles.axes1,data.freqString);


totalTones = length(get(handles.listbox2,'String'));

cla(handles.axes2);
for t = 1:totalTones
    response = data.toneResponses{neuron,attn}(t,:);
    if t == tone
        width = 4;
    else
        width = 0.5;
    end
    plot(handles.axes2,response,'LineWidth',width);
    xticks(handles.axes2,0:30:180);
    xticklabels(handles.axes2,{'-1','0','1','2','3','4','5'});
    xlabel(handles.axes2,'Time after stimulus onset (sec)');
    ylabel(handles.axes2,'\Delta F/std(F)');
    hold on;
end

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


% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2
global data;

neuron = get(handles.listbox1,'Value');
attn = get(handles.listbox3,'Value');
tone = get(handles.listbox2,'Value');

cla(handles.axes1);
for a = 1:length(data.attenuations)
    tCurve = squeeze(data.tuningCurves(neuron,a,:));

    semilogx(handles.axes1,data.frequencies,tCurve);
    xlabel(handles.axes1,'Frequency (Hz)');
    ylabel(handles.axes1,'\Delta F/std(F)');
    hold(handles.axes1,'on');
end
legend(handles.axes1,data.attnString);

% ylim(handles.axes1,[-0.2 1]);
% set(handles.axes1,'xlim',[3 5],'xticks',[3 4 5],'xticklabels',{1 2 3},'xlabel','Frequency (kHz)');
% xticks(handles.axes1,[3 4 5]
% xticklabels(handles.axes1,data.freqString);


totalTones = length(get(handles.listbox2,'String'));

cla(handles.axes2);
for t = 1:totalTones
    response = data.toneResponses{neuron,attn}(t,:);
    if t == tone
        width = 4;
    else
        width = 0.5;
    end
    plot(handles.axes2,response,'LineWidth',width);
    xticks(handles.axes2,0:30:180);
    xticklabels(handles.axes2,{'-1','0','1','2','3','4','5'});
    xlabel(handles.axes2,'Time after stimulus onset (sec)');
    ylabel(handles.axes2,'\Delta F/std(F)');
    hold on;
end


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
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

[file folder] = uigetfile('/Users/Aaron/Documents');
data = load([folder '/' file]);
for f = 1:length(data.frequencies)
    freq{f} = num2str(data.frequencies(f));
end
data.freqString = freq;

set(handles.text4,'String',file);

for i = 1:size(data.tuningCurves,1)
    list1{i} = num2str(i);
end
set(handles.listbox1,'String',list1);

for j = 1:size(data.tuningCurves,3)
    list2{j} = num2str(j);
end
set(handles.listbox2,'String',list2);

for a = 1:length(data.attenuations)
    attn{a} = num2str(70 + data.attenuations(a));
end
data.attnString = attn;
set(handles.listbox3,'String',attn);


% --- Executes on selection change in listbox3.
function listbox3_Callback(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox3
global data;

neuron = get(handles.listbox1,'Value');
attn = get(handles.listbox3,'Value');
tone = get(handles.listbox2,'Value');

cla(handles.axes1);
for a = 1:length(data.attenuations)
    tCurve = squeeze(data.tuningCurves(neuron,a,:));

    semilogx(handles.axes1,data.frequencies,tCurve);
    xlabel(handles.axes1,'Frequency (Hz)');
    ylabel(handles.axes1,'\Delta F/std(F)');
    hold(handles.axes1,'on');
end
legend(handles.axes1,data.attnString);

% ylim(handles.axes1,[-0.2 1]);
% set(handles.axes1,'xlim',[3 5],'xticks',[3 4 5],'xticklabels',{1 2 3},'xlabel','Frequency (kHz)');
% xticks(handles.axes1,[3 4 5]
% xticklabels(handles.axes1,data.freqString);


totalTones = length(get(handles.listbox2,'String'));

cla(handles.axes2);
for t = 1:totalTones
    response = data.toneResponses{neuron,attn}(t,:);
    if t == tone
        width = 4;
    else
        width = 0.5;
    end
    plot(handles.axes2,response,'LineWidth',width);
    xticks(handles.axes2,0:30:180);
    xticklabels(handles.axes2,{'-1','0','1','2','3','4','5'});
    xlabel(handles.axes2,'Time after stimulus onset (sec)');
    ylabel(handles.axes2,'\Delta F/std(F)');
    hold on;
end


% --- Executes during object creation, after setting all properties.
function listbox3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end