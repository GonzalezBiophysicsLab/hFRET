function varargout = yaxisprompt(varargin)
% yaxisprompt helps the user construct a tree for hhmm_gui_main
% the output here gets shunted directly into hhmm_gui_main
% JH 9/21/2016

% Last Modified by GUIDE v2.5 28-Mar-2018 09:32:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @yaxisprompt_OpeningFcn, ...
                   'gui_OutputFcn',  @yaxisprompt_OutputFcn, ...
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


% --- Executes just before yaxisprompt is made visible.
function yaxisprompt_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to yaxisprompt (see VARARGIN)

% Choose default command line output for yaxisprompt
handles.input = varargin{1};
if handles.input ~= 0
    set(handles.promptext,'string','Set the max limit');
end
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes yaxisprompt wait for user response (see UIRESUME)
handles.outputs.bottom = 'a';
handles.outputs.top = 'b';
guidata(hObject, handles);
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = yaxisprompt_OutputFcn(hObject, eventdata, handles) 
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
% try
    b = double(str2double(a));
    if ~isfloat(handles.outputs.bottom)&&handles.input==0
        handles.outputs.bottom = b;
        handles.outputs.count = 1;
        if isfloat(handles.outputs.bottom)
            set(handles.promptext,'string','Set the max limit');
        end
    elseif ~isfloat(handles.outputs.top)
        handles.outputs.top = b;
        handles.outputs.count = 1;
        if isfloat(handles.outputs.top)            
            uiresume(handles.figure1);
            set(handles.promptext,'string','generating....');
        end
    end
    guidata(hObject, handles);
% end
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
