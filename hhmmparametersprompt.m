function varargout = hhmmparametersprompt(varargin)
% hhmmparametersprompt helps the user construct a tree for hhmm_gui_main
% the output here gets shunted directly into hhmm_gui_main
% JH 9/21/2016

% Last Modified by GUIDE v2.5 28-Mar-2018 10:23:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @hhmmparametersprompt_OpeningFcn, ...
                   'gui_OutputFcn',  @hhmmparametersprompt_OutputFcn, ...
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


% --- Executes just before hhmmparametersprompt is made visible.
function hhmmparametersprompt_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to hhmmparametersprompt (see VARARGIN)

% Choose default command line output for hhmmparametersprompt
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes hhmmparametersprompt wait for user response (see UIRESUME)
handles.outputs.maxIter = 'a';
handles.outputs.restarts = 'b';
handles.outputs.baseline = 'c';
handles.outputs.guess = 'd';
handles.outputs.ampfrac = 'e';
guidata(hObject, handles);
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = hhmmparametersprompt_OutputFcn(hObject, eventdata, handles) 
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
    handles.outputs.response = a;
    if ~isfloat(handles.outputs.maxIter)
        handles.outputs.maxIter = b;
        handles.outputs.count = 1;
        if isfloat(handles.outputs.maxIter)
            set(handles.promptext,'string','Set the maximum number of restarts');
        end
    elseif ~isfloat(handles.outputs.restarts)
        handles.outputs.restarts = b;
        if isfloat(handles.outputs.restarts)
            set(handles.promptext,'string','Shared emissions? (y/n)');
        end
    elseif (strcmp(handles.outputs.response,'y')||strcmp(handles.outputs.response,'Y'))&&~islogical(handles.outputs.baseline)
        handles.outputs.baseline = 1==0;
        set(handles.promptext,'string','Automatic emissions? (y/n)');
    elseif (strcmp(handles.outputs.response,'n')||strcmp(handles.outputs.response,'N'))&&~islogical(handles.outputs.baseline)
        handles.outputs.baseline = 1==1;
        set(handles.promptext,'string','Automatic emissions? (y/n)');
    elseif (strcmp(handles.outputs.response,'y')||strcmp(handles.outputs.response,'Y'))&&strcmp(handles.outputs.guess,'d')
        handles.outputs.guess = 'auto';
        set(handles.promptext,'string','Analyze fractional amplitudes? (y/n)');
    elseif (strcmp(handles.outputs.response,'n')||strcmp(handles.outputs.response,'N'))&&strcmp(handles.outputs.guess,'d')
        set(handles.promptext,'string','Enter guess (e.g.: [1 2 3 4]). Number of guesses must equal number of centers.');
    elseif strcmp(handles.outputs.guess,'d')
        tmp = sort(str2num(a));
        handles.outputs.guess = tmp(:);
        set(handles.promptext,'string','Analyze fractional amplitudes? (y/n)');
    elseif (strcmp(handles.outputs.response,'y')||strcmp(handles.outputs.response,'Y'))&&strcmp(handles.outputs.ampfrac,'e')
        handles.outputs.ampfrac = 0;
        uiresume(handles.figure1);
        set(handles.promptext,'string','generating....');
    elseif (strcmp(handles.outputs.response,'n')||strcmp(handles.outputs.response,'N'))&&strcmp(handles.outputs.ampfrac,'e')
        handles.outputs.ampfrac = 0;
        set(handles.promptext,'string','Specify amplitude channel (e.g.: 1)');
    elseif handles.outputs.ampfrac == 0
        handles.outputs.ampfrac = b;
        uiresume(handles.figure1);
        set(handles.promptext,'string','generating....');
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
