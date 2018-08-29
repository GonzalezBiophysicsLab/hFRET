function varargout = vbscope_sort_traces(varargin)
% sorts traces whose output matches vbscope_analyze_traces formats
% 
% Users must manually input their sort preferences into the table


% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help vbscope_sort_traces

% Last Modified by GUIDE v2.5 06-Sep-2016 09:38:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @vbscope_sort_traces_OpeningFcn, ...
                   'gui_OutputFcn',  @vbscope_sort_traces_OutputFcn, ...
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


% --- Executes just before vbscope_sort_traces is made visible.
function vbscope_sort_traces_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to vbscope_sort_traces (see VARARGIN)

% Choose default command line output for vbscope_sort_traces
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes vbscope_sort_traces wait for user response (see UIRESUME)


if ~isempty(varargin)
    grow = varargin{1};
    handles.tracestats = grow.tracestats;
    handles.class_logical = ones(10,1);
    if isfield(grow.tracestats,'order')
        for rk = 1:10
            handles.class_string{rk,1} = num2str(grow.tracestats.order(rk,:));
        end
    else
        handles.class_string = {'1';'2';'3';'4';'5';'6';'7';'8';'9';'10'};
    end
    for rk = 0:1:9
        handles.class_string{rk+1,2} = num2str(sum(handles.tracestats.sort==rk));
    end 
    set(handles.classshow,'Data',handles.class_string,'ColumnEditable',logical([1 0]));
end
handles.sorttype_logical = 1;
guidata(hObject, handles);

uiwait(handles.figure1);




% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);

% --- Outputs from this function are returned to the command line.
function varargout = vbscope_sort_traces_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.tracestats;
delete(hObject);

% --- Executes on button press in applyordering.
function applyordering_Callback(hObject, eventdata, handles)
% hObject    handle to applyordering (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.sorttype_logical == 1
    handles.class_string = get(handles.classshow,'data');
    for rk = 1:1:10
        order(rk,:) = str2num(handles.class_string{rk,1});
    end
    % get the raw ordering
    [index,loc] = unique(order);
    class_index = 0:1:9;
    % apply to class indices
    grow = handles.tracestats;
    grow = rmfield(grow,'data');
    grow = rmfield(grow,'sort');
    grow = rmfield(grow,'manual');
    grow.order = order;
    cn = grow.colornum;
    % prepare a struct to replace the old tracestats
    for zk = 1:10
        class = handles.tracestats.sort==class_index(loc(zk));
        % yes if in class; no if not
        if zk == 1
            grow.sort = handles.tracestats.sort(class);
            grow.manual = handles.tracestats.manual(class);
            for rrk = 1:cn
                grow.data.colors(rrk).traces = handles.tracestats.data.colors(rrk).traces(:,class);
                grow.data.colors(rrk).spots = handles.tracestats.data.colors(rrk).spots(:,class);
                grow.data.colors(rrk).excluded = handles.tracestats.data.colors(1).excluded(class,:);
                grow.data.colors(rrk).start = handles.tracestats.data.colors(rrk).start(class,:);
                grow.data.colors(rrk).photobleach = handles.tracestats.data.colors(rrk).photobleach(class,:);
            end
            % populate the first
        else
            grow.sort = [grow.sort;handles.tracestats.sort(class)];
            grow.manual = [grow.manual;handles.tracestats.manual(class)];
            for rrk = 1:cn
                grow.data.colors(rrk).traces = [grow.data.colors(rrk).traces handles.tracestats.data.colors(rrk).traces(:,class)];
                grow.data.colors(rrk).spots = [grow.data.colors(rrk).spots handles.tracestats.data.colors(rrk).spots(:,class)];
                grow.data.colors(rrk).excluded = [grow.data.colors(1).excluded;handles.tracestats.data.colors(1).excluded(class,:)];
                grow.data.colors(rrk).start = [grow.data.colors(rrk).start;handles.tracestats.data.colors(rrk).start(class,:)];
                grow.data.colors(rrk).photobleach = [grow.data.colors(rrk).photobleach;handles.tracestats.data.colors(rrk).photobleach(class,:)];
            end
            % grow the subsequent
        end
    end
    handles.tracestats = grow;
elseif handles.sorttype_logical == 2
    cn = handles.tracestats.colornum;
    corr = zeros(size(handles.tracestats.data.colors(1).traces,2),1);
    for rk = 1:cn
        for rrk = rk:cn
            if rk~=rrk
                for zk = 1:size(handles.tracestats.data.colors(1).traces,2)
                    trace1 = diff(handles.tracestats.data.colors(rk).traces(:,zk));
                    trace2 = diff(handles.tracestats.data.colors(rrk).traces(:,zk));
                    temp = ifft(fft(trace1).*conj(fft(trace2)));
                    corr(zk,:) = corr(zk,:)+temp(1);
                end
            end
        end
    end
    [~,rule] = sort(corr);
    handles.tracestats.sort = handles.tracestats.sort(rule);
    handles.tracestats.manual = handles.tracestats.manual(rule);
    for rrk = 1:cn
        handles.tracestats.data.colors(rrk).traces = handles.tracestats.data.colors(rrk).traces(:,rule);
        handles.tracestats.data.colors(rrk).spots = handles.tracestats.data.colors(rrk).spots(:,rule);
        handles.tracestats.data.colors(1).excluded = handles.tracestats.data.colors(1).excluded(rule,:);
        handles.tracestats.data.colors(rrk).start = handles.tracestats.data.colors(rrk).start(rule,:);
        handles.tracestats.data.colors(rrk).photobleach = handles.tracestats.data.colors(rrk).photobleach(rule,:);     
    end
end
guidata(hObject, handles);
uiresume(handles.figure1);
% hand back and close the gui


% --- Executes on key press with focus on classshow and none of its controls.
function classshow_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to classshow (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in sorttype.
function sorttype_Callback(hObject, eventdata, handles)
% hObject    handle to sorttype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns sorttype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from sorttype


stri = get(hObject,'String');
val = get(hObject,'Value');

switch stri{val};
    case 'Sort by class          '
        handles.sorttype_logical = 1;
    case 'Sort by anticorrelation'
        handles.sorttype_logical = 2;
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function sorttype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sorttype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
