function varargout = coordinateprompt(varargin)
% coordinateprompt helps the user construct a tree for hhmm_gui_main
% the output here gets shunted directly into hhmm_gui_main
% JH 9/21/2016

% Last Modified by GUIDE v2.5 04-Apr-2018 13:41:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @coordinateprompt_OpeningFcn, ...
                   'gui_OutputFcn',  @coordinateprompt_OutputFcn, ...
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


% --- Executes just before coordinateprompt is made visible.
function coordinateprompt_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to coordinateprompt (see VARARGIN)

% Choose default command line output for coordinateprompt
handles.depth = varargin{1};
handles.instructions = varargin{2};
if strcmp(handles.instructions,'coordplot')
    set(handles.promptext,'string','Direct only? (y/n)');
else
    
end
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes coordinateprompt wait for user response (see UIRESUME)
guidata(hObject, handles);

if numel(handles.depth) == 1
    handles.outputs = 'Direct';
    guidata(hObject, handles);
else
    uiwait(handles.figure1);
end


% --- Outputs from this function are returned to the command line.
function varargout = coordinateprompt_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.outputs;
delete(hObject);



function input_Callback(hObject, eventdata, handles)
% hObject    handle to input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input as text
%        str2double(get(hObject,'String')) returns contents of input as a double

a = get(hObject,'String');
if strcmp(a,'y')||strcmp(a,'Y')&&strcmp(handles.instructions,'coordplot')
    handles.instructions = 'done';
    handles.outputs = 'Direct';
    set(handles.promptext,'string','generating....');
elseif strcmp(a,'n')||strcmp(a,'N')&&strcmp(handles.instructions,'coordplot')
    set(handles.promptext,'string','Set dimension depth');
    handles.instructions = 1;
elseif handles.instructions == 1
    set(handles.promptext,'string','generating....');
    handles.instructions = 'done';
    handles.outputs = str2double(a);
end
guidata(hObject, handles);

if strcmp(handles.instructions,'done')
    uiresume(handles.figure1);
end
% --- Executes during object creation, after setting all properties.
function input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
