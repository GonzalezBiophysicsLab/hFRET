function varargout = port(varargin)
% PORT MATLAB code for port.fig
%      PORT, by itself, creates a new PORT or raises the existing
%      singleton*.
%
%      H = PORT returns the handle to a new PORT or the handle to
%      the existing singleton*.
%
%      PORT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PORT.M with the given input arguments.
%
%      PORT('Property','Value',...) creates a new PORT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before port_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to port_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help port

% Last Modified by GUIDE v2.5 27-Aug-2018 19:51:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @port_OpeningFcn, ...
                   'gui_OutputFcn',  @port_OutputFcn, ...
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


% --- Executes just before port is made visible.
function port_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to port (see VARARGIN)

% Choose default command line output for port
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes port wait for user response (see UIRESUME)
% uiwait(handles.figure1);

if ~isempty(varargin)
    grow = varargin{1};
    handles.tracestats = grow.tracestats;
    handles.plotstuff = grow.plotstuff;
    handles.class_string = {'On';'On';'On';'On';'On';'On';'On';'On';'On';'On'};
    handles.class_logical = ones(10,1);
    for rk = 0:1:9
        handles.class_string{rk+1,2} = num2str(sum(handles.tracestats.sort==rk));
    end
    set(handles.classshow,'Data',handles.class_string);
end
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = port_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in analysisprogram.
function analysisprogram_Callback(hObject, eventdata, handles)
% hObject    handle to analysisprogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns analysisprogram contents as cell array
%        contents{get(hObject,'Value')} returns selected item from analysisprogram

stri = get(hObject,'String');
val = get(hObject,'Value');

switch stri{val};
    case 'vbFRET'
        set(handles.instructions,'string',{'(1) Select a save location for the .dat file containing the classes you want to analyze';'';'(2) Select the path to the upper level of vbFRET (folder containing src)'});        
        cn = handles.tracestats.colornum;
        
        logi = [handles.class_logical(end);handles.class_logical(1:end-1)];
        for rk = cn:-1:1
            temp = [];
            for zk = 1:length(logi)
                if logi(zk)==1
                    class = handles.tracestats.sort==(zk-1);
                    temp = [temp handles.tracestats.data.colors(rk).traces(:,class)];
                end
            end
            t = 1:cn:size(temp,2)*cn;
            data(:,t+rk-1) = temp;
        end
        [fn,pn] = uiputfile('*.dat');
        disp('Saving')
        dlmwrite([pn fn],data);
        disp('Saved')
        path = uigetdir;
        addpath(genpath(path));
        vbFRET;
    case 'ebFRET'
        set(handles.instructions,'string',{'(1) Select a save location for the .dat file containing the classes you want to analyze';'';'(2) Select the path to the upper level of ebFRET (folder containing src)'});        
        cn = handles.tracestats.colornum;
        
        logi = [handles.class_logical(end);handles.class_logical(1:end-1)];
        for rk = cn:-1:1
            temp = [];
            for zk = 1:length(logi)
                if logi(zk)==1
                    class = handles.tracestats.sort==(zk-1);
                    temp = [temp handles.tracestats.data.colors(rk).traces(:,class)];
                end
            end
            t = 1:cn:size(temp,2)*cn;
            data(:,t+rk-1) = temp;
        end
        [fn,pn] = uiputfile('*.dat');
        disp('Saving')
        dlmwrite([pn fn],data);
        disp('Saved')
        path = uigetdir;
        addpath(genpath([path '/src']));
        ebf = ebFRET();
    case 'BIASD'
    case 'Subpopulation HMM'
    case 'HFRET'
        cn = handles.tracestats.colornum;
        logi = [handles.class_logical(end);handles.class_logical(1:end-1)];
        class = handles.tracestats.sort~=handles.tracestats.sort;
        for zk = 1:length(logi)
            if logi(zk)==1
                class = class|handles.tracestats.sort==(zk-1);
            end
        end
            handles.tracestats.sort = handles.tracestats.sort(class);
            handles.tracestats.manual = handles.tracestats.manual(class);
        for rk = cn:-1:1
            handles.tracestats.data.colors(rk).traces = handles.tracestats.data.colors(rk).traces(:,class);
            handles.tracestats.data.colors(rk).excluded = handles.tracestats.data.colors(rk).excluded(class);
            handles.tracestats.data.colors(rk).start = handles.tracestats.data.colors(rk).start(class);
            handles.tracestats.data.colors(rk).photobleach = handles.tracestats.data.colors(rk).photobleach(class);
        end
        handles.tracestats.tracenum = sum(class);
        out.tracestats = handles.tracestats;
        out.plotstuff = handles.plotstuff;
        hfret_gui_main(out);
end
% --- Executes during object creation, after setting all properties.
function analysisprogram_CreateFcn(hObject, eventdata, handles)
% hObject    handle to analysisprogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on classshow and none of its controls.
function classshow_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to classshow (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)



s = eventdata.Key;

for rk = 1:10
    if rk == 10
        rrk = 0;
        class = 1;
    else
        rrk = rk;
        class = rk+1;
    end
    if strcmp(s,num2str(rrk))
        handles.class_logical(rk,:) = ~handles.class_logical(rk,:);
        if handles.class_logical(rk,:)
            handles.class_string{class,1} = 'On';
        else
            handles.class_string{class,1} = 'Off';            
        end
        guidata(hObject,handles);
    end
end
set(handles.classshow,'Data',handles.class_string);
