function varargout = treeprompt(varargin)
% treeprompt helps the user construct a tree for hhmm_gui_main
% the output here gets shunted directly into hhmm_gui_main
% JH 9/21/2016

% Last Modified by GUIDE v2.5 21-Sep-2016 13:31:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @treeprompt_OpeningFcn, ...
                   'gui_OutputFcn',  @treeprompt_OutputFcn, ...
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


% --- Executes just before treeprompt is made visible.
function treeprompt_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to treeprompt (see VARARGIN)

% Choose default command line output for treeprompt
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes treeprompt wait for user response (see UIRESUME)
handles.outputs.depth = 'c';
guidata(hObject, handles);
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = treeprompt_OutputFcn(hObject, eventdata, handles) 
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
try
    b = uint8(str2double(a));
    if ~isinteger(handles.outputs.depth)
        handles.outputs.depth = b;
        handles.outputs.count = 1;
        if isinteger(handles.outputs.depth)
            set(handles.promptext,'string',['Set the number of values for level ' num2str(handles.outputs.count) '/' num2str(handles.outputs.depth)]);
        end
    else
        if handles.outputs.depth==1&&handles.outputs.count>=2
            handles.outputs.k = double(b);
            handles.outputs.count = handles.outputs.count+1;
        else
            handles.outputs.depth_vector(handles.outputs.count) = b;
            handles.outputs.count = handles.outputs.count+1;            
        end
        if handles.outputs.count>handles.outputs.depth
            if handles.outputs.depth ~= 1
                uiresume(handles.figure1);
                set(handles.promptext,'string','generating....');
            elseif handles.outputs.depth == 1&&handles.outputs.count>2
                uiresume(handles.figure1);
                set(handles.promptext,'string','generating....');
            else
                set(handles.promptext,'string','Set the number of static subpopulations');
            end

        elseif isinteger(handles.outputs.depth)
            set(handles.promptext,'string',['Set the number of values for level ' num2str(handles.outputs.count) '/' num2str(handles.outputs.depth)]);
        end
    end
    guidata(hObject, handles);
catch
    set(handles.promptext,'string',['Error, try again; Set the number of values for level ' num2str(handles.outputs.count) '/' num2str(handles.outputs.depth)]);
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
