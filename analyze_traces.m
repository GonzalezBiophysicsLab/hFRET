function varargout = analyze_traces(varargin)
% ANALYZE_TRACES MATLAB code for analyze_traces.fig
%      ANALYZE_TRACES, by itself, creates a new ANALYZE_TRACES or raises the existing
%      singleton*.
%
%      H = ANALYZE_TRACES returns the handle to a new ANALYZE_TRACES or the handle to
%      the existing singleton*.
%
%      ANALYZE_TRACES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ANALYZE_TRACES.M with the given input arguments.
%
%      ANALYZE_TRACES('Property','Value',...) creates a new ANALYZE_TRACES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before analyze_traces_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to analyze_traces_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help analyze_traces

% Last Modified by GUIDE v2.5 27-Aug-2018 19:51:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @analyze_traces_OpeningFcn, ...
                   'gui_OutputFcn',  @analyze_traces_OutputFcn, ...
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




% --- Executes just before analyze_traces is made visible.
function analyze_traces_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to analyze_traces (see VARARGIN)

% Choose default command line output for analyze_traces
handles.output = hObject;

handles.figures.spots(1).object = handles.color1;
handles.figures.spots(2).object = handles.color2;
handles.figures.spots(3).object = handles.color3;
handles.figures.spots(4).object = handles.color4;
for rk = 1:length(handles.figures.spots)
    axes(handles.figures.spots(rk).object);
    axis off;
end
axes(handles.histogram);
axis off;


if ~isempty(varargin)
    grow = varargin{1};
    handles.tracestats.colornum = length(grow.colors);
    cn = handles.tracestats.colornum;
    
    if cn == 1
        handles.tracestats.colornum = 1;
        handles.plotstuff.colors = [0 .5 0];
    elseif cn == 2
        handles.tracestats.colornum = 2;
        handles.plotstuff.colors = [0 .5 0;1 0 0];
    elseif cn == 3
        handles.tracestats.colornum = 3;
        handles.plotstuff.colors = [0 .5 0;1 0 0;0 0 .5];
    elseif cn == 4
        handles.tracestats.colornum = 4;
        handles.plotstuff.colors = [0 .5 0;1 0 0;0 0 .5;0 0 0];
    end
    for rrk = 1:cn
        if ~isfield(grow.colors(rrk),'excluded')
            grow.colors(rrk).excluded = zeros(size(grow.colors(rrk).traces,2),1);
        end
    end
    for rk = cn:-1:1
        fun1 = @(x) length(grow.colors(rk).traces(grow.colors(rk).traces(:)>x))-round(.99*numel(grow.colors(rk).traces));
        fun2 = @(x) length(grow.colors(rk).traces(grow.colors(rk).traces(:)<x))-round(.99*numel(grow.colors(rk).traces));
        mm(rk,:) = fzero(fun1,mean(grow.colors(rk).traces(:)));
        m(rk,:) = fzero(fun2,mean(grow.colors(rk).traces(:)));
    end

    handles.plotstuff.yaxis = [mean(mm)/2 mean(m)]*6;
    handles.plotstuff.xaxis = [1 length(grow.colors(1).spots(:,1))];
    handles.plotstuff.saveaxes = [handles.plotstuff.yaxis handles.plotstuff.xaxis];
    handles.plotstuff.caxis = [mean(mm)/2 mean(m)]*6;
    handles.plotstuff.current_trace = 1;
    handles.plotstuff.current_time = 1;
    handles.tracestats.tracenum = length(grow.colors(1).excluded);
    handles.tracestats.data = grow;
    for rk = 1:cn
        handles.tracestats.data.colors(rk).start = ones(handles.tracestats.tracenum,1);
        handles.tracestats.data.colors(rk).photobleach = size(handles.tracestats.data.colors(rk).traces,1)*ones(handles.tracestats.tracenum,1);
    end
    handles.tracestats.sort = zeros(handles.tracestats.tracenum,1);
    handles.plotstuff.plotlim = round(mean(m)/1e2)*1e2;
    handles.plotstuff.photobleaching = 1==0;
    handles.tracestats.threshold = 0.95;
    handles.tracestats.threshold_up = 0.95;
    handles.tracestats.manual = zeros(handles.tracestats.tracenum,1);
    guidata(hObject,handles);

    x = handles.plotstuff.xaxis(end);
    set(handles.exploretraces,'sliderstep',[1/(handles.tracestats.tracenum-1),10/(handles.tracestats.tracenum-1)],'max',1,'min',1/(handles.tracestats.tracenum),'Value',1/(handles.tracestats.tracenum));
    set(handles.adjustdot,'sliderstep',[1/(x-1),10/(x-1)],'max',1,'min',1/x,'Value',1/x);
    set(handles.brightness,'sliderstep',[.01,.25],'max',2,'min',0,'Value',1/2);
    set(handles.pbslider,'sliderstep',[.0005 .0005],'max',1,'min',0.5,'Value',0.95);
    set(handles.pbslider2,'sliderstep',[.0005 .0005],'max',1,'min',0.5,'Value',0.95);
    scopeplot(handles);
    disp('Data has been loaded')
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes analyze_traces wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = analyze_traces_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in analyze.
function analyze_Callback(hObject, eventdata, handles)
% hObject    handle to analyze (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns analyze contents as cell array
%        contents{get(hObject,'Value')} returns selected item from analyze


% --- Executes during object creation, after setting all properties.
function analyze_CreateFcn(hObject, eventdata, handles)
% hObject    handle to analyze (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in exclude.
function exclude_Callback(hObject, eventdata, handles)
% hObject    handle to exclude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of exclude
check = get(hObject,'Value');
a = handles.plotstuff.current_trace;
cn = handles.tracestats.colornum;
for rk = 1:cn
    if check
        handles.tracestats.data.colors(rk).excluded(a,:) = 1;
    end
end
guidata(hObject,handles);

% --- Executes on slider movement.
function adjustdot_Callback(hObject, eventdata, handles)
% hObject    handle to adjustdot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


x = round(get(hObject,'Value')*handles.plotstuff.xaxis(end));
cn = handles.tracestats.colornum;
a = handles.plotstuff.current_trace;
handles.plotstuff.current_time = x;
guidata(hObject,handles);
for rk = 1:cn
    axes(handles.figures.spots(rk).object);
    if ~isempty(handles.tracestats.data.colors(rk).spots{x,a})
        attempt = 0;
        ok = 0;
        while ~ok
            try
                imagesc(double(reshape(handles.tracestats.data.colors(rk).spots{x,a},5,5-attempt)));
                ok = 1;
            catch
                attempt = attempt+1;
            end
            if attempt>5
                ok = 1;
            end
        end
        colormap('hot')
        caxis(handles.plotstuff.caxis)
    end
    axis off
end
axes(handles.intensity)
cla
hold on
for rk = 1:cn
    plot(handles.tracestats.data.colors(rk).traces(:,a),'color',handles.plotstuff.colors(rk,:))
end
plot([x x],handles.plotstuff.yaxis,'k--')
xlim(handles.plotstuff.xaxis)
ylim(handles.plotstuff.yaxis)
hold off
axes(handles.fret)
cla
maxint = zeros(size(handles.tracestats.data.colors(1).traces,1),1);
hold on
if cn > 1
    for rk = 1:cn
        maxint = maxint+handles.tracestats.data.colors(rk).traces(:,a);
    end
    for rk = 2:cn
        plot(handles.tracestats.data.colors(rk).traces(:,a)./maxint,'color',handles.plotstuff.colors(rk,:));
    end
    ylim([-.75 2])
    xlim([1 length(maxint)])
    plot([x x],[-.75 2],'k--')
else
    plot([x x],[-.75 2],'k--')
    axis off
end
hold off

% --- Executes during object creation, after setting all properties.
function adjustdot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to adjustdot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in loaddata.
function loaddata_Callback(hObject, eventdata, handles)
% hObject    handle to loaddata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename,pathname,~] = uigetfile('*','MultiSelect','on');
cn = handles.tracestats.colornum;
disp('Loading data');
if iscell(filename)
    for rk = 1:length(filename)
        z = load([pathname filename{rk}]);
        if isstruct(z)
            if rk == 1
                grow = z.vbscope_outputs;
                for rrk = 1:cn
                    if ~isfield(grow.colors(rrk),'spots')
                        grow.colors(rrk).spots = cell(size(grow.colors(rrk).traces,1),size(grow.colors(rrk).traces,2));
                    end
                    if ~isfield(grow.colors(rrk),'excluded')
                        grow.colors(rrk).excluded = zeros(size(grow.colors(rrk).traces,2),1);
                    end
                end
            else
                for rrk = 1:cn
                    grow.colors(rrk).traces = [grow.colors(rrk).traces z.vbscope_outputs.colors(rrk).traces];
                    if isfield(z.vbscope_outputs.colors(rrk),'spots')
                        grow.colors(rrk).spots = [grow.colors(rrk).spots z.vbscope_outputs.colors(rrk).spots];
                    else
                        grow.colors(rrk).spots = [grow.colors(rrk).spots cell(size(z.vbscope_outputs.colors(rrk).traces,1),size(z.vbscope_outputs.colors(rrk).traces,2))];
                    end
                    if isfield(z.vbscope_outputs.colors(rrk),'excluded')
                        grow.colors(rrk).excluded = [grow.colors(rrk).excluded;z.vbscope_outputs.colors(rrk).excluded];
                    else
                        grow.colors(rrk).excluded = [grow.colors(rrk).excluded;zeros(size(z.vbscope_outputs.colors(rrk).traces,2),1)];
                    end
                end
            end
        elseif ismatrix(z)
            if rk == 1
                t = 1:cn:size(z,2);
                for rrk = 1:cn
                    grow.colors(rrk).traces = z(:,t+(rrk-1));
                    grow.colors(rrk).spots = cell(size(z,1),size(z,2)/cn);
                    grow.colors(rrk).excluded = zeros(size(z,2)/cn,1);
                end
            else
                t = 1:cn:size(z,2);
                for rrk = 1:cn
                    grow.colors(rrk).traces = [grow.colors(rrk).traces z(:,t+(rrk-1))];
                    grow.colors(rrk).spots = [grow.colors(rrk).spots cell(size(z,1),size(z,2)/cn)];
                    grow.colors(rrk).excluded = [grow.colors(rrk).excluded;zeros(size(z,2)/cn,1)];
                end
            end
        else
            disp('Illegal input');
        end
    end
else
    z = load([pathname filename]);
    if isstruct(z)
        grow = z.vbscope_outputs;
        if ~isfield(grow.colors(1),'spots')
            for rrk = 1:cn
                grow.colors(rrk).spots = cell(size(grow.colors(rrk).traces,1),size(grow.colors(rrk).traces,2));
            end
        end
        for rrk = 1:cn
            if ~isfield(grow.colors(rrk),'excluded')
                grow.colors(rrk).excluded = zeros(size(grow.colors(rrk).traces,2),1);
            end
        end
    else
        t = 1:cn:size(z,2);
        for rrk = 1:cn
            grow.colors(rrk).traces = z(:,t+(rrk-1));
            grow.colors(rrk).spots = cell(size(z,1),size(z,2)/cn);
            grow.colors(rrk).excluded = zeros(size(z,2)/cn,1);
        end
    end
end
for rk = cn:-1:1
    fun1 = @(x) length(grow.colors(rk).traces(grow.colors(rk).traces(:)>x))-round(.99*numel(grow.colors(rk).traces));
    fun2 = @(x) length(grow.colors(rk).traces(grow.colors(rk).traces(:)<x))-round(.99*numel(grow.colors(rk).traces));
    mm(rk,:) = fzero(fun1,mean(grow.colors(rk).traces(:)));
    m(rk,:) = fzero(fun2,mean(grow.colors(rk).traces(:)));
end

handles.plotstuff.yaxis = [mean(mm)/2 mean(m)]*6;
handles.plotstuff.xaxis = [1 length(grow.colors(1).spots(:,1))];
handles.plotstuff.saveaxes = [handles.plotstuff.yaxis handles.plotstuff.xaxis];
handles.plotstuff.caxis = [mean(mm)/2 mean(m)]*6;
handles.plotstuff.current_trace = 1;
handles.plotstuff.current_time = 1;
if isfield(handles.tracestats,'data')
    for rrk = 1:cn
        handles.tracestats.data.colors(rrk).traces = [handles.tracestats.data.colors(rrk).traces grow.colors(rrk).traces];
        handles.tracestats.data.colors(rrk).spots = [handles.tracestats.data.colors(rrk).spots grow.colors(rrk).spots];
        handles.tracestats.data.colors(rrk).excluded = [handles.tracestats.data.colors(rrk).excluded;grow.colors(rrk).excluded];
    end
else
    handles.tracestats.data = grow;
end
handles.tracestats.tracenum = length(handles.tracestats.data.colors(1).excluded);
for rk = 1:cn
    handles.tracestats.data.colors(rk).start = ones(handles.tracestats.tracenum,1);
    handles.tracestats.data.colors(rk).photobleach = size(handles.tracestats.data.colors(rk).traces,1)*ones(handles.tracestats.tracenum,1);
end
handles.tracestats.sort = zeros(handles.tracestats.tracenum,1);
handles.plotstuff.plotlim = round(mean(m)/1e2)*1e2;
handles.plotstuff.photobleaching = 1==0;
handles.tracestats.threshold = 0.95;
handles.tracestats.threshold_up = 0.95;
handles.tracestats.manual = zeros(handles.tracestats.tracenum,1);
guidata(hObject,handles);

x = handles.plotstuff.xaxis(end);
set(handles.exploretraces,'sliderstep',[1/(handles.tracestats.tracenum-1),10/(handles.tracestats.tracenum-1)],'max',1,'min',1/(handles.tracestats.tracenum),'Value',1/(handles.tracestats.tracenum));
set(handles.adjustdot,'sliderstep',[1/(x-1),10/(x-1)],'max',1,'min',1/x,'Value',1/x);
set(handles.brightness,'sliderstep',[.01,.25],'max',1,'min',0,'Value',1/2);
set(handles.pbslider,'sliderstep',[.0005 .0005],'max',1,'min',0.5,'Value',0.95);
set(handles.pbslider2,'sliderstep',[.0005 .0005],'max',1,'min',0.5,'Value',0.95);
scopeplot(handles);
disp('Data has been loaded')




% --- Executes on button press in reanalyze_dot.
function reanalyze_dot_Callback(hObject, eventdata, handles)
% hObject    handle to reanalyze_dot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in colornum_popup.
function colornum_popup_Callback(hObject, eventdata, handles)
% hObject    handle to colornum_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns colornum_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from colornum_popup



stri = get(hObject,'String');
val = get(hObject,'Value');



switch stri{val};
    case '1 color'
        handles.tracestats.colornum = 1;
        handles.plotstuff.colors = [0 .5 0];
    case '2 color'
        handles.tracestats.colornum = 2;
        handles.plotstuff.colors = [0 .5 0;1 0 0];
    case '3 color'
        handles.tracestats.colornum = 3;
        handles.plotstuff.colors = [0 .5 0;1 0 0;0 0 .5];
    case '4 color'
        handles.tracestats.colornum = 4;
        handles.plotstuff.colors = [0 .5 0;1 0 0;0 0 .5;0 0 0];
end
guidata(hObject,handles);




% --- Executes during object creation, after setting all properties.
function colornum_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to colornum_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function exploretraces_Callback(hObject, eventdata, handles)
% hObject    handle to exploretraces (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider



a = round(get(hObject,'Value')*handles.tracestats.tracenum);
handles.plotstuff.current_trace = a;
set(handles.sortclass,'string',['Class ' num2str(handles.tracestats.sort(a))]);
guidata(hObject,handles);
scopeplot(handles);
% --- Executes during object creation, after setting all properties.
function exploretraces_CreateFcn(hObject, eventdata, handles)
% hObject    handle to exploretraces (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function brightness_Callback(hObject, eventdata, handles)
% hObject    handle to brightness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

a = get(hObject,'Value')*handles.plotstuff.plotlim*2;
handles.plotstuff.yaxis(end) = a;
handles.plotstuff.caxis(end) = a/4;
cn = handles.tracestats.colornum;

axes(handles.intensity)
ylim(handles.plotstuff.yaxis);

for rk = 1:cn
    axes(handles.figures.spots(rk).object);
    caxis(handles.plotstuff.caxis);
end
guidata(hObject,handles);






% --- Executes during object creation, after setting all properties.
function brightness_CreateFcn(hObject, eventdata, handles)
% hObject    handle to brightness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in savetraces.
function savetraces_Callback(hObject, eventdata, handles)
% hObject    handle to savetraces (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cn = handles.tracestats.colornum;
t = 1:cn:(handles.tracestats.tracenum-sum(handles.tracestats.data.colors(1).excluded))*cn;
if handles.savestuff.type == 1
    [fn,pn] = uiputfile('*.dat');
    handles.tracestats.data.colors(1).excluded = zeros(size(handles.tracestats.data.colors(1).start,1),1);
    for rk = cn:-1:1
        out(:,t+rk-1) = handles.tracestats.data.colors(rk).traces(:,~handles.tracestats.data.colors(1).excluded);
    end
    disp('Saving')
    dlmwrite([pn fn],out);
    disp('Saved')
else
    [fn,pn] = uiputfile('*.mat');
    handles.tracestats.data.colors(1).excluded = zeros(size(handles.tracestats.data.colors(1).start,1),1);
    for rk = cn:-1:1
        vbscope_outputs.colors(rk).traces = handles.tracestats.data.colors(rk).traces;
        vbscope_outputs.colors(rk).spots = handles.tracestats.data.colors(rk).spots;
        vbscope_outputs.colors(rk).excluded = handles.tracestats.data.colors(1).excluded;
    end
    disp('Saving')
    save([pn fn],'-v7.3','vbscope_outputs');
    disp('Saved')
end

% --- Executes on selection change in filetype.
function filetype_Callback(hObject, eventdata, handles)
% hObject    handle to filetype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns filetype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from filetype



stri = get(hObject,'String');
val = get(hObject,'Value');

switch stri{val};
    case 'dat type file'
        handles.savestuff.type = 1;
    case 'vbscope type file'
        handles.savestuff.type = 2;
end
guidata(hObject,handles);




% --- Executes during object creation, after setting all properties.
function filetype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filetype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in dotevolution.
function dotevolution_Callback(hObject, eventdata, handles)
% hObject    handle to dotevolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

a = handles.plotstuff.current_trace;
cn = handles.tracestats.colornum;
for zk = 1:handles.plotstuff.xaxis(end)
    for rk = 1:cn
        axes(handles.figures.spots(rk).object);
        cla
        if ~isempty(handles.tracestats.data.colors(rk).spots{zk,a})
            attempt = 0;
            ok = 0;
            while ~ok
                try
                    imagesc(double(reshape(handles.tracestats.data.colors(rk).spots{zk,a},5,5-attempt)));
                    ok = 1;
                catch
                    attempt = attempt+1;
                end
                if attempt>5
                    ok = 1;
                end
            end
            caxis(handles.plotstuff.caxis)
            colormap('hot')
        end
        axis off
    end
    pause(0.001);
end


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

s = eventdata.Key;

for rk = 0:1:9
    if strcmp(s,num2str(rk))
        handles.tracestats.sort(handles.plotstuff.current_trace) = rk;
        set(handles.sortclass,'string',['Class ' num2str(rk)]);
        guidata(hObject,handles);
        disp(['Class ' num2str(rk)])
    end
end

switch s
    case 'n'
        set(handles.exploretraces,'Value',max([1/handles.tracestats.tracenum min([1 (handles.plotstuff.current_trace+1)/handles.tracestats.tracenum])]));

        handles.plotstuff.current_time = 1;

        a = round(max([1/handles.tracestats.tracenum min([1 (handles.plotstuff.current_trace+1)/handles.tracestats.tracenum])])*handles.tracestats.tracenum);
        handles.plotstuff.current_trace = a;
        set(handles.sortclass,'string',['Class ' num2str(handles.tracestats.sort(a))]);
        guidata(hObject,handles);
        scopeplot(handles);
    case 'b'
        set(handles.exploretraces,'Value',max([1/handles.tracestats.tracenum min([1 (handles.plotstuff.current_trace-1)/handles.tracestats.tracenum])]));
        handles.plotstuff.current_time = 1;

        a = round(max([1/handles.tracestats.tracenum min([1 (handles.plotstuff.current_trace-1)/handles.tracestats.tracenum])])*handles.tracestats.tracenum);
        handles.plotstuff.current_trace = a;
        set(handles.sortclass,'string',['Class ' num2str(handles.tracestats.sort(a))]);
        guidata(hObject,handles);
        scopeplot(handles);
    case 'f'
        
        cn = handles.tracestats.colornum;
        tl = size(handles.tracestats.data.colors(1).traces,1);
        a = handles.plotstuff.current_trace;
        handles.tracestats.manual(a) = 1;
        for rk = 1:cn
            [x,~] = getpts(handles.intensity);
            if isempty(x)||x<0||x>tl
                x = tl;
            end
            handles.tracestats.data.colors(rk).photobleach(a) = round(x);
        end
        guidata(hObject,handles);
        scopeplot(handles);
    case 's'
        
        cn = handles.tracestats.colornum;
        tl = size(handles.tracestats.data.colors(1).traces,1);
        a = handles.plotstuff.current_trace;
        for rk = 1:cn
            [x,~] = getpts(handles.intensity);
            if isempty(x)||x<0||x>tl
                x = 1;
            end
            handles.tracestats.data.colors(rk).start(a) = round(x);
        end
        guidata(hObject,handles);
        scopeplot(handles);
    case 'z'
        zoom;

end


% --- Executes on button press in sortag.
function sortag_Callback(hObject, eventdata, handles)
% hObject    handle to sortag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename,pathname,~] = uiputfile('*.dat');
handles.tracestats.data.colors(1).excluded = zeros(size(handles.tracestats.data.colors(1).start,1),1);
e = ~handles.tracestats.data.colors(1).excluded;
output = [];
for rk = 1:handles.tracestats.colornum
    output = [output handles.tracestats.data.colors(rk).start(e) handles.tracestats.data.colors(rk).photobleach(e)];
end
dlmwrite([pathname filename],[handles.tracestats.sort(e) output]);


% --- Executes on button press in estpb.
function estpb_Callback(hObject, eventdata, handles)
% hObject    handle to estpb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.plotstuff.photobleaching = ~handles.plotstuff.photobleaching;
guidata(hObject,handles);
cn = handles.tracestats.colornum;
tn = handles.tracestats.tracenum;
tl = size(handles.tracestats.data.colors(1).traces,1);
factor = 10;
t = 1:factor:tl;
if handles.plotstuff.photobleaching
    for zk = 1:cn
        k = [];
        eff = [];
        for rk = 1:tn
            for rrk = 1:factor
                temp2 = handles.tracestats.data.colors(zk).traces(rrk:end,rk);
                temp2(temp2<0) = 0;
                a = diff(downsample(temp2,factor));
                b = t*0;
                b(1:length(a)) = a;
                k(t+rrk-1,rk) = b;
            end
            eff = [eff;handles.tracestats.data.colors(zk).traces(:,rk)];
        end
        [~,p] = kmeans(eff,2);
        m = reshape(k,numel(k),1);
        best = scope_vmp_gauss_mix_mod(m/1e3,3,[-abs(diff(p)) 0 abs(diff(p))]*3/4/1e3);
        down = reshape(best.conprob(:,1),size(k,1),size(k,2));
        up = reshape(best.conprob(:,3),size(k,1),size(k,2));
        for rk = tn:-1:1
            if ~handles.tracestats.manual(rk)
                b = find(down(:,rk)>handles.tracestats.threshold,1,'last');
                a = find(up(:,rk)>handles.tracestats.threshold_up,1,'last');
                if isempty(a)
                    a = 0;
                end
                if isempty(b)||b<1.5||a>b
                    temp(rk,:) = tl;
                else
                    temp(rk,:) = b-1;
                end
            else
                temp(rk,:) = handles.tracestats.data.colors(zk).photobleach(rk,:);
            end
        end
        handles.tracestats.data.colors(zk).photobleach = temp;
        handles.tracestats.data.pb(zk).probabilities.down = down;
        handles.tracestats.data.pb(zk).probabilities.up = up;
    end
    guidata(hObject,handles);
    scopeplot(handles);
end
% Hint: get(hObject,'Value') returns toggle state of estpb


% --- Executes on button press in manualfix.
function manualfix_Callback(hObject, eventdata, handles)
% hObject    handle to manualfix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cn = handles.tracestats.colornum;
tl = size(handles.tracestats.data.colors(1).traces,1);
a = handles.plotstuff.current_trace;
handles.tracestats.manual(a) = 1;
for rk = 1:cn
    [x,~] = getpts(handles.intensity);
    if isempty(x)||x<0||x>tl
        x = tl;
    end
    handles.tracestats.data.colors(rk).photobleach(a) = round(x);
end
guidata(hObject,handles);

scopeplot(handles);


% --- Executes on button press in manualstart.
function manualstart_Callback(hObject, eventdata, handles)
% hObject    handle to manualstart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


cn = handles.tracestats.colornum;
tl = size(handles.tracestats.data.colors(1).traces,1);
a = handles.plotstuff.current_trace;
for rk = 1:cn
    [x,~] = getpts(handles.intensity);
    if isempty(x)||x<0||x>tl
        x = 1;
    end
    handles.tracestats.data.colors(rk).start(a) = round(x);
end
guidata(hObject,handles);

scopeplot(handles);


function out = scopeplot(handles)

out = [];
a = handles.plotstuff.current_trace;
set(handles.tracenumtext,'string',['Trace ' num2str(a) '/' num2str(handles.tracestats.tracenum)]);

cn = handles.tracestats.colornum;
handles.plotstuff.current_time = 1;
tl = size(handles.tracestats.data.colors(1).traces,1);

if handles.plotstuff.photobleaching
    for rk = cn:-1:1
        pbtime(:,rk) = handles.tracestats.data.colors(rk).photobleach(a,:);
    end
else
    pbtime = ones(1,cn)*(tl-1);
end
pbtime(pbtime>tl-1) = tl-1;
for rk = cn:-1:1
    start(rk,:) = handles.tracestats.data.colors(rk).start(a,:);
end
axes(handles.intensity)
cla
hold on
for rk = 1:cn
    plot(start(rk):pbtime(rk)+1,handles.tracestats.data.colors(rk).traces(start(rk):pbtime(rk)+1,a),'color',handles.plotstuff.colors(rk,:))
    plot(pbtime(rk):tl,handles.tracestats.data.colors(rk).traces(pbtime(rk):end,a),':','MarkerSize',10,'color',handles.plotstuff.colors(rk,:))
    plot(1:start(rk),handles.tracestats.data.colors(rk).traces(1:start(rk),a),':','MarkerSize',10,'color',handles.plotstuff.colors(rk,:))
end
xlim(handles.plotstuff.xaxis)
ylim(handles.plotstuff.yaxis)
hold off
axes(handles.fret)
cla
maxint = zeros(size(handles.tracestats.data.colors(1).traces,1),1);
hold on
if cn > 1
    for rk = 1:cn
        maxint = maxint+handles.tracestats.data.colors(rk).traces(:,a);
    end
    for rk = 2:cn
        fret(:,rk) = handles.tracestats.data.colors(rk).traces(:,a)./maxint;
        plot(max(start):min(pbtime)+1,fret(max(start):min(pbtime)+1,rk),'color',handles.plotstuff.colors(rk,:));
        plot(min(pbtime)-1:tl,fret(min(pbtime)-1:end,rk),':','MarkerSize',10,'color',handles.plotstuff.colors(rk,:));
        plot(1:max(start),fret(1:max(start),rk),':','MarkerSize',10,'color',handles.plotstuff.colors(rk,:));
    end
    ylim([-.75 2])
    xlim(handles.plotstuff.xaxis)
else
    axis off
end
hold off

hold on
for rk = 1:cn
    axes(handles.figures.spots(rk).object);
    if ~isempty(handles.tracestats.data.colors(rk).spots{1,a})
        attempt = 0;
        ok = 0;
        while ~ok
            try
                imagesc(double(reshape(handles.tracestats.data.colors(rk).spots{1,a},5,5-attempt)));
                ok = 1;
            catch
                attempt = attempt+1;
            end
            if attempt>5
                ok = 1;
            end
        end
        colormap('hot')
        caxis(handles.plotstuff.caxis)
    end
    axis off
end
hold off
axes(handles.histogram)

hold on
cla


for rk = 2:cn
    temp = fret(max(start):min(pbtime),rk);
    xx = -.75:(2.75)/(4*length(temp)^.33):2;
    nn = xx*0;
    for zk = 1:length(xx)-1
        nn(zk) = sum(temp>xx(zk)&temp<xx(zk+1));
    end
    plot(handles.histogram,nn/max(nn),xx,'color',handles.plotstuff.colors(rk,:))
    ylim([-.75 2])
    axis off
end
hold off


% --- Executes on selection change in analysismenu.
function analysismenu_Callback(hObject, eventdata, handles)
% hObject    handle to analysismenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns analysismenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from analysismenu



stri = get(hObject,'String');
val = get(hObject,'Value');
output.tracestats = handles.tracestats;
output.plotstuff = handles.plotstuff;

switch stri{val};
    case 'Model-Free'
        vbscope_nonparam_analysis(output);
    case 'Parametric'
        port(output);
end
guidata(hObject,handles);





% --- Executes during object creation, after setting all properties.
function analysismenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to analysismenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function pbslider_Callback(hObject, eventdata, handles)
% hObject    handle to pbslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

a = get(hObject,'Value');
handles.tracestats.threshold = a;
set(handles.pbtext,'string',['Downward transition cutoff: ' num2str(a)]);
tn = handles.tracestats.tracenum;
tl = size(handles.tracestats.data.colors(1).traces,1);

try
    for zk = 1:handles.tracestats.colornum
        down = handles.tracestats.data.pb(zk).probabilities.down;
        up = handles.tracestats.data.pb(zk).probabilities.up;
        for rk = tn:-1:1
            if ~handles.tracestats.manual
                b = find(down(:,rk)>handles.tracestats.threshold,1,'last');
                a = find(up(:,rk)>handles.tracestats.threshold_up,1,'last');
                if isempty(a)
                    a = 0;
                end
                if isempty(b)||b<1.5||a>b
                    temp(rk,:) = tl;
                else
                    temp(rk,:) = b-1;
                end
            end
        end
        handles.tracestats.data.colors(zk).photobleach = temp;
    end
end
scopeplot(handles);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function pbslider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pbslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function pbslider2_Callback(hObject, eventdata, handles)
% hObject    handle to pbslider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

a = get(hObject,'Value');
handles.tracestats.threshold_up = a;
set(handles.pbtext2,'string',['Upward transition cutoff: ' num2str(a)]);
tn = handles.tracestats.tracenum;
tl = size(handles.tracestats.data.colors(1).traces,1);

try
    for zk = 1:handles.tracestats.colornum
        down = handles.tracestats.data.pb(zk).probabilities.down;
        up = handles.tracestats.data.pb(zk).probabilities.up;
        for rk = tn:-1:1
            b = find(down(:,rk)>handles.tracestats.threshold,1,'last');
            a = find(up(:,rk)>handles.tracestats.threshold_up,1,'last');
            if isempty(a)
                a = 0;
            end
            if isempty(b)||b<1.5||a>b
                temp(rk,:) = tl;
            else
                temp(rk,:) = b-1;
            end
        end
        handles.tracestats.data.colors(zk).photobleach = temp;
    end
end
scopeplot(handles);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function pbslider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pbslider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in loadindices.
function loadindices_Callback(hObject, eventdata, handles)
% hObject    handle to loadindices (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename,pathname,~] = uigetfile('*','MultiSelect','on');
if iscell(filename)
    get = [];
    for rk = 1:length(filename)
        z = load([pathname filename{rk}]);
        get = [get;z];
    end
else
    get = load([pathname filename]);
end

handles.tracestats.sort = get(:,1);
count = 1;
for rk = 2:2:handles.tracestats.colornum*2
    handles.tracestats.data.colors(count).start = get(:,rk);
    handles.tracestats.data.colors(count).photobleach = get(:,rk+1);
    count = count+1;
end
handles.plotstuff.photobleaching = ~handles.plotstuff.photobleaching;
set(handles.estpb,'Value',1);
guidata(hObject,handles);
scopeplot(handles);


% --- Executes on button press in zoomstop.
function zoomstop_Callback(hObject, eventdata, handles)
% hObject    handle to zoomstop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

zoom off
handles.plotstuff.yaxis = ylim(handles.intensity);
handles.plotstuff.xaxis = xlim(handles.intensity);
handles.plotstuff.xaxis(1) = max([1 handles.plotstuff.xaxis(1)]);
scopeplot(handles);
guidata(hObject,handles);


% --- Executes on button press in resetzoom.
function resetzoom_Callback(hObject, eventdata, handles)
% hObject    handle to resetzoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

zoom off
handles.plotstuff.yaxis = handles.plotstuff.saveaxes(1:2);
handles.plotstuff.xaxis = handles.plotstuff.saveaxes(3:4);
guidata(hObject,handles);
scopeplot(handles);


% --- Executes on button press in sorttraces.
function sorttraces_Callback(hObject, eventdata, handles)
% hObject    handle to sorttraces (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.tracestats = vbscope_sort_traces(handles);
guidata(hObject,handles);
scopeplot(handles);

a = handles.plotstuff.current_trace;
set(handles.sortclass,'string',['Class ' num2str(handles.tracestats.sort(a))]);
