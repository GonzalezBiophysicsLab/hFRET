function varargout = hfret_gui_main(varargin)
% includes trace splash, 2 popups for different kinds of outputs, error
% analysis, and some ensemble level plots
%

% Last Modified by GUIDE v2.5 27-Aug-2018 19:50:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @hfret_gui_main_OpeningFcn, ...
                   'gui_OutputFcn',  @hfret_gui_main_OutputFcn, ...
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


% --- Executes just before hfret_gui_main is made visible.
function hfret_gui_main_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to hfret_gui_main (see VARARGIN)

% Choose default command line output for hfret_gui_main
handles.output = hObject;

if ~isempty(varargin)
    handles.tracestats = varargin{1}.tracestats;
    handles.plotstuff = varargin{1}.plotstuff;
end
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes hfret_gui_main wait for user response (see UIRESUME)
% uiwait(handles.graph1);

axes(handles.graph1);
axis off;
axes(handles.graph2);
axis off;
axes(handles.intensity);
axis off;
axes(handles.coordinates);
axis off;
axes(handles.fret);
axis off;
axes(handles.histogram);
axis off;
axes(handles.modelscheme);
axis off;
warning('off','all');
handles.traceposition = 1;
set(handles.traceslider,'sliderstep',[1/(handles.tracestats.tracenum-1),10/(handles.tracestats.tracenum-1)],'max',1,'min',1/(handles.tracestats.tracenum),'Value',1/(handles.tracestats.tracenum));
set(handles.modelslider,'sliderstep',[1/(10-1),10/(10-1)],'max',1,'min',1/(10),'Value',1/(10));
try
    set(handles.xaxisslider,'sliderstep',[1/(10-1),10/(10-1)],'max',1,'min',1/(10),'Value',0);
end
set(handles.timetext,'string',[]);
% this has an invalid setting so that it will not be rendered.

handles.hhmm_current = 1;
for rk = 1:10
    handles.hhmm{rk}.evidence = -Inf;
end
set(handles.modeltext,'string',['Model ' num2str(1)]);
handles.plotstuff.coordlog = 1==0;
handles.plotstuff.current_trace = 1;
handles.plotstuff.hxaxis = handles.plotstuff.xaxis;
scopeplot(handles);
guidata(hObject,handles);

% --- Outputs from this function are returned to the command line.
function varargout = hfret_gui_main_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function traceslider_Callback(hObject, eventdata, handles)
% hObject    handle to traceslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


a = round(get(hObject,'Value')*handles.tracestats.tracenum);
handles.plotstuff.current_trace = a;
set(handles.class,'string',['Class ' num2str(handles.tracestats.sort(a))]);
guidata(hObject,handles);
scopeplot(handles);

% --- Executes during object creation, after setting all properties.
function traceslider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to traceslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on selection change in analysistype.
function analysistype_Callback(hObject, eventdata, handles)
% hObject    handle to analysistype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns analysistype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from analysistype


% --- Executes during object creation, after setting all properties.
function analysistype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to analysistype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in savetype.
function savetype_Callback(hObject, eventdata, handles)
% hObject    handle to savetype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns savetype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from savetype


% --- Executes during object creation, after setting all properties.
function savetype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to savetype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function hhmm_plot(handles)
    handles.traceposition = a;


% --- Executes on selection change in plottype.
function plottype_Callback(hObject, eventdata, handles)
% hObject    handle to plottype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns plottype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plottype


% --- Executes during object creation, after setting all properties.
function plottype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plottype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function population_Callback(hObject, eventdata, handles)
% hObject    handle to population (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function population_CreateFcn(hObject, eventdata, handles)
% hObject    handle to population (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in graphdesign.
function graphdesign_Callback(hObject, eventdata, handles)
% hObject    handle to graphdesign (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

out = treeprompt;
handles.hhmm_options.depth_vector = out.depth_vector;
try
    handles.hhmm_options.k = out.k;
catch
    handles.hhmm_options.k = 1;
end
plotdots(handles.modelscheme,handles.hhmm_options.depth_vector,handles.hhmm_options.k);
handles.hhmm_options.production = double(handles.hhmm_options.depth_vector(end));
handles.hhmm_options.depth_vector = double(handles.hhmm_options.depth_vector(1:end-1));
guidata(hObject,handles);

function plotdots(hax,depth,k)
    axes(hax)
    cla(hax)
    hold on
    depth = [1;depth(:)];
    for nk = 0:k-1
        t = (1:2:(prod(depth))*2);
        for rk = 2:length(depth)
            rrk = length(depth)-rk+1;
            z = mean(t(1:(prod(depth(1:rrk+1))/prod(depth(1:rrk)))))/2;
            tt = 2*z*(1:2:(2*prod(depth(1:rrk))));
            count = 1;
            if rk ~= length(depth)
                plot(hax,tt+nk*max(tt),0*tt+rk,'.k','MarkerSize',20)
            end
            for zzk = 1:length(t)
                plot(hax,[t(zzk)+nk*max(t) tt(count)+nk*max(t)],[rk-1 rk],'-k')
                if mod(zzk,depth(rrk+1))==0
                    count = count+1; 
                end
            end
            t = tt;
        end
        t = 1:2:(prod(depth))*2;
        plot(hax,t+nk*max(t),0*(t)+1,'.r','MarkerSize',20)
    end    
    hold off
    ylim(hax,[0 length(depth)+1])
    xlim(hax,[0 (prod(depth))*2+(k-1)*(prod(depth))*2])
    axis off

function colordots(hax,colors,depth,k,input,logij)
    axes(hax)
    hold on
    depth = [1;depth(:)];
    t = (1:2:(prod(depth))*2);
    rk = length(depth)-input;
    
    for nk = 0:k-1
        if input < length(depth)-1
            if rk>2
                for rrrk = 2:rk-1
                    rrk = length(depth)-rrrk+1;
                    z = mean(t(1:(prod(depth(1:rrk+1))/prod(depth(1:rrk)))))/2;
                    tt = 2*z*(1:2:(2*prod(depth(1:rrk))));
                    t = tt;
                end
            end
            rrk = length(depth)-rk+1;
            z = mean(t(1:(prod(depth(1:rrk+1))/prod(depth(1:rrk)))))/2;
            tt = 2*z*(1:2:(2*prod(depth(1:rrk))));
            for zk = 1:length(tt)
                plot(hax,tt(zk)+nk*max(tt),rk,'.','MarkerSize',20,'color',colors(zk,:))
            end
        else
            for zk = 1:length(t)
                if logij
                    zzk = depth(end)-mod(zk,depth(end));
                else
                    zzk = zk;
                end
                plot(hax,t(zk)+nk*max(t),1,'.','MarkerSize',20,'color',colors(zzk,:))
            end
        end
    end    
    hold off
    ylim(hax,[0 length(depth)+1])
    xlim(hax,[0 (prod(depth))*2+(k-1)*(prod(depth))*2])
    axis off
   
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
if cn~=1
    axes(handles.intensity)
    axis on
else
    axes(handles.intensity)
    axis off
    axes(handles.fret)
    axis on
end
ylabel('Amplitude')
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
if cn~=1
    axes(handles.fret)
    axis on
    ylabel('Fractional Amplitude')
    cla
    maxint = zeros(size(handles.tracestats.data.colors(1).traces,1),1);
    hold on
    if cn > 1
        for rk = 1:cn
            maxint = maxint+handles.tracestats.data.colors(rk).traces(:,a);
        end
        for rk = 2:cn
            fret(:,rk) = handles.tracestats.data.colors(rk).traces(:,a)./maxint;
            plot(start(1):min(pbtime)+1,fret(start(1):min(pbtime)+1,rk),'color',handles.plotstuff.colors(rk,:));
            plot(min(pbtime)-1:tl,fret(min(pbtime)-1:end,rk),':','MarkerSize',10,'color',handles.plotstuff.colors(rk,:));
            plot(1:start(1),fret(1:start(1),rk),':','MarkerSize',10,'color',handles.plotstuff.colors(rk,:));
        end
        ylim([-.75 2])
        xlim(handles.plotstuff.xaxis)
    else
        axis off
    end
    hold off
end
axes(handles.histogram)
hold on
cla

if cn==1
    rk = 1;
    temp = handles.tracestats.data.colors(rk).traces(start(rk):pbtime(rk)+1,a);
    xx = handles.plotstuff.yaxis(1):2.75/(4*length(temp)^.33):handles.plotstuff.yaxis(2);
    nn = xx*0;
    for zk = 1:length(xx)-1
        nn(zk) = sum(temp>xx(zk)&temp<xx(zk+1));
    end
    plot(handles.histogram,nn/max(nn),xx,'color',handles.plotstuff.colors(rk,:))
    ylim(handles.plotstuff.yaxis)
else   
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
end

axis off
hold off

axes(handles.coordinates)
if handles.plotstuff.coordlog
    b = handles.hhmm_current;
    axis on
    cla
    if strcmp(handles.plotstuff.coordlevel,'Direct')
        time = start(1):min(pbtime)+1;
        hold on
        if isempty(handles.hhmm{b}.inputoptions.depth_vector)
            temp = handles.hhmm{b}.vbem{a}.ideals;
        else
            temp = handles.hhmm{b}.ideals{a};
        end        
        if handles.hhmm{b}.inputoptions.ampfrac~=0
            plot(handles.coordinates,time,temp,'-k','LineWidth',.5);
        else
            plot(handles.coordinates,time,temp/1e2,'-k','LineWidth',.5);
        end
        for rk = 1:handles.hhmm{b}.inputoptions.production
            if handles.hhmm{b}.inputoptions.baseline
                try
                    mu = handles.hhmm{b}.emissions(a).m;
                catch
                    mu = handles.hhmm{b}.vbem{a}.m;
                end
            else
                try
                    mu = handles.hhmm{b}.emissions.m;
                catch
                    mu = handles.hhmm{b}.m;
                end
            end
            map = temp==mu(rk);
            if handles.hhmm{b}.inputoptions.ampfrac~=0
                plot(handles.coordinates,time(map),temp(map),'.','MarkerSize',9,'color',handles.plotstuff.coordcolors(rk,:))
            else
                plot(handles.coordinates,time(map),temp(map)/1e2,'.','MarkerSize',9,'color',handles.plotstuff.coordcolors(rk,:))
            end
        end
        hold off
        if handles.hhmm{b}.inputoptions.ampfrac~=0
            ylabel('Amplitude')
            ylim(handles.plotstuff.yaxis)
        else
            ylabel('Fractional Amplitude')
            ylim([-.75 2])
        end
        set(handles.coordinates,'yticklabelmode','auto');
    else
        ylabel('Coordinate')
        time = start(1):min(pbtime)+1;
        temp = handles.hhmm{b}.pathstate{a}(:,handles.plotstuff.coordlevel+1);
        hold on
        plot(handles.coordinates,time,temp,'-k','LineWidth',1);
        for rk = 1:size(handles.plotstuff.coordcolors,1)
            map = temp==rk;
            plot(handles.coordinates,time(map),temp(map),'.','MarkerSize',10,'color',handles.plotstuff.coordcolors(rk,:))
        end
        hold off
        set(handles.coordinates,'yticklabelmode','auto');
        ylim([0,prod(handles.hhmm{b}.depth_vector(1:handles.plotstuff.coordlevel))+1])
    end
    xlim(handles.plotstuff.xaxis)
else
    axis off
    cla
end


% --- Executes on button press in analyze.
function analyze_Callback(hObject, eventdata, handles)
% hObject    handle to analyze (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cn = handles.tracestats.colornum;
last = handles.tracestats.tracenum;
data = cell(last,1);
handles.plotstuff.current_time = 1;
tl = size(handles.tracestats.data.colors(1).traces,1);

if cn == 1
    count = 1;
    for rrk = 1:last
        for rk = cn:-1:1
            pbtime(:,rk) = handles.tracestats.data.colors(rk).photobleach(rrk,:);
        end
        pbtime(pbtime>tl-1) = tl-1;
        for rk = cn:-1:1
            start(rk,:) = handles.tracestats.data.colors(rk).start(rrk,:);
        end
        data{count} = handles.tracestats.data.colors(1).traces(start(1):min(pbtime)+1,rrk);
        count = count+1;
        if handles.hhmm_options.baseline&&strcmp(handles.hhmm_options.guess,'auto')
            [~,m] = kmeans(data{count-1},handles.hhmm_options.production);
            tmp(:,count-1) = sort(m);
        elseif handles.hhmm_options.baseline&&~strcmp(handles.hhmm_options.guess,'auto')
            tmp(:,count-1) = sort(handles.hhmm_options.guess);
        end
    end
end
% load up data into acceptable format


if cn~=1
    if handles.hhmm_options.ampfrac == 0
        maxint = zeros(size(handles.tracestats.data.colors(1).traces,1),1);
        for rrk = 1:last
            for rk = cn:-1:1
                pbtime(:,rk) = handles.tracestats.data.colors(rk).photobleach(rrk,:);
            end
            pbtime(pbtime>tl-1) = tl-1;
            for rk = cn:-1:1
                start(rk,:) = handles.tracestats.data.colors(rk).start(rrk,:);
            end
            for rk = 1:cn
                maxint = maxint+handles.tracestats.data.colors(rk).traces(:,rrk);
            end
            for rk = 2:cn
                fret(:,rk) = handles.tracestats.data.colors(rk).traces(:,rrk)./maxint;
            end
            data{rrk} = fret(start(1):min(pbtime)+1,2)*1e2;
            % scale
            if handles.hhmm_options.baseline&&strcmp(handles.hhmm_options.guess,'auto')
                [~,m] = kmeans(data{rrk},handles.hhmm_options.production);
                tmp(:,rrk) = sort(m);
            elseif handles.hhmm_options.baseline&&~strcmp(handles.hhmm_options.guess,'auto')
                tmp(:,rrk) = sort(handles.hhmm_options.guess);
            end
        end
    else
        for rrk = 1:last
            for rk = cn:-1:1
                pbtime(:,rk) = handles.tracestats.data.colors(rk).photobleach(rrk,:);
            end
            pbtime(pbtime>tl-1) = tl-1;
            for rk = cn:-1:1
                start(rk,:) = handles.tracestats.data.colors(rk).start(rrk,:);
            end
            data{rrk} = handles.tracestats.data.colors(handles.hhmm_options.ampfrac).traces(start(1):min(pbtime)+1,rrk);
            if handles.hhmm_options.baseline&&strcmp(handles.hhmm_options.guess,'auto')
                [~,m] = kmeans(data{rrk},handles.hhmm_options.production);
                tmp(:,rrk) = sort(m);
            elseif handles.hhmm_options.baseline&&~strcmp(handles.hhmm_options.guess,'auto')
                tmp(:,rrk) = sort(handles.hhmm_options.guess);
            end
        end
    end
end
if (~handles.hhmm_options.baseline)&&strcmp(handles.hhmm_options.guess,'auto')
    [~,m] = kmeans(cell2mat(data(:)),handles.hhmm_options.production);
    tmp = sort(m);
elseif ~handles.hhmm_options.baseline
    tmp = handles.hhmm_options.guess;
end
if handles.hhmm_options.ampfrac == 0
    if handles.hhmm_options.baseline
        tmp = tmp*1e2;
    else
        tmp = handles.hhmm_options.guess(:)*1e2;
    end
end
% get the options

if isempty(handles.hhmm_options.depth_vector)
    rate = zeros(length(data),handles.hhmm_options.production*(handles.hhmm_options.production-1));
    vb = cell(length(data),1);
    disp('Getting prior');
    for rk = 1:length(data)
        if handles.hhmm_options.baseline
            vb{rk} = vmp_hmm_modular(data{rk},handles.hhmm_options.production,tmp(:,rk));
        else
            vb{rk} = vmp_hmm_modular(data{rk},handles.hhmm_options.production,handles.hhmm_options.guess(:));
        end
        try
            temp = vb{rk}.mark.mix';
            rate(rk,:) = temp(~eye(handles.hhmm_options.production,handles.hhmm_options.production));
        end
    end
    handles.hhmm{handles.hhmm_current} = vmp_xi_consensus_smfret(data,vb,rate,handles.hhmm_options.k,handles.hhmm_options.production,handles.hhmm_options.maxIter,handles.hhmm_options.baseline,tmp);
    handles.hhmm{handles.hhmm_current}.depth_vector = handles.hhmm_options.production;
else
    options = handles.hhmm_options;
    options.guess = tmp;
    handles.hhmm{handles.hhmm_current} = vmp_hhmm(data,options);
    handles.hhmm{handles.hhmm_current}.depth_vector = [handles.hhmm_options.depth_vector handles.hhmm_options.production];
end
handles.hhmm{handles.hhmm_current}.inputoptions = handles.hhmm_options;
handles.hhmm{handles.hhmm_current}.data = data;
handles.plotstuff.coordlog = 1==0;
guidata(hObject,handles);
disp('Analysis complete')

function dwellData = hhmm_getdwells(paths)
% standard in the Gonzalez lab; get dwell times

dwellData = [];

for p=(1:length(paths))
    try
        path = paths{p};
        % transition happens on frames when next state is different from current state
        % and on last frame
        transitions = find([path(1:end-1) ~= path(2:end); true]); 
        lengths = transitions - [0; transitions(1:end-1)];
        fret = path(transitions);
        nextFret = [path(transitions(1:end-1)+1); NaN];
        dwellData = [dwellData; p*ones(size(fret)), fret, nextFret, lengths];
    end
end

function [map,t] = fxn_hist2D(scaled,paths,a,b,limits)

count = 1;
while ~isempty(a)
    try
        a = a(a(:,2)~=0,:);
        b = b(b(:,2)~=0,:);
        if ~isempty(a)
            spot = find((b(:,2)-a(1,2))>0&b(:,1)==a(1,1),1);
            ok = 1;
            if isempty(spot)
                spot = find(b(:,1)~=a(1,1),1);
                ok = 0;
            else
                try
                    prev = a(1,2)-find(paths{a(1,1)}(1:(a(1,2)-2))~=paths{a(1,1)}(a(1,2)-1),1,'last');
                end
                prev = min([prev limits(3)]);
                sync{count} = scaled{a(1,1)}((a(1,2)-prev:b(spot,2)));
                len2(count,:) = b(spot,2)-a(1);
                if ~isempty(sync{count})
                    count = count+1;
                end
            end
            if ok == 1
                a(a(:,2)<b(spot,2)&a(1,1)==a(:,1),2) = 0;
                b(a(:,2)<b(spot,2)&a(1,1)==a(:,1),2) = 0;
            else
                a = a(2:end,:);
            end
        end
    catch
        a = a(2:end,:);
    end
end
sync = sync(1:end-1);

for rk = length(sync):-1:1
    len(rk,:) = length(sync{rk});
end
sync_together = zeros(max(len),length(sync));
for rk = 1:length(sync)
    sync_together(1:length(sync{rk}),rk) = sync{rk};
end

t = limits(1):abs(diff(limits))/1e2:limits(2);
map = zeros(length(t),max([size(sync_together,1) 1e3]));

for zk = 1:length(t)-1
    map(zk,1:size(sync_together,1)) = sum(sync_together>t(zk)&sync_together<t(zk+1),2);
end

% --- Executes on button press in yaxis.
function yaxis_Callback(hObject, eventdata, handles)
% hObject    handle to yaxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


out = yaxisprompt(0);
handles.plotstuff.yaxis = [out.bottom out.top];
scopeplot(handles);
guidata(hObject,handles);

% --- Executes on button press in hhmmparameters.
function hhmmparameters_Callback(hObject, eventdata, handles)
% hObject    handle to hhmmparameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


out = hhmmparametersprompt;
handles.hhmm_options.maxIter = out.maxIter;
handles.hhmm_options.restarts = out.restarts;
handles.hhmm_options.baseline = out.baseline;
handles.hhmm_options.guess = out.guess;
handles.hhmm_options.ampfrac = out.ampfrac;
guidata(hObject,handles);


% --- Executes on slider movement.
function modelslider_Callback(hObject, eventdata, handles)
% hObject    handle to modelslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

a = round(get(hObject,'Value')*10);
handles.hhmm_current = a;
set(handles.modeltext,'string',['Model ' num2str(a)]);
try
    handles.hhmm_options = handles.hhmm{a}.inputoptions;
    plotdots(handles.modelscheme,handles.hhmm{a}.depth_vector,handles.hhmm_options.k);
    handles.plotstuff.coordlog = 1==0;
catch
    handles.modelscheme;
    handles.plotstuff.coordlog = 1==0;
    cla
end
scopeplot(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function modelslider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to modelslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in savebutton.
function savebutton_Callback(hObject, eventdata, handles)
% hObject    handle to savebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fn,pn] = uiputfile('*.mat');
n = handles.hhmm_current;
out = handles.hhmm{n};
save([pn fn],'-v7.3','out');
disp('Saved!');

% --- Executes on button press in comparebutton.
function comparebutton_Callback(hObject, eventdata, handles)
% hObject    handle to comparebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

evidence = zeros(10,1);
depths = cell(10,1);
for rk = 1:10
    try
        try
            evidence(rk,:) = handles.hhmm{rk}.record.evidence(end);
        catch
            evidence(rk,:) = handles.hhmm{rk}.record(end);
        end
        depths{rk} = handles.hhmm{rk}.depth_vector;
        k{rk} = handles.hhmm{rk}.inputoptions.k;
    catch
        evidence(rk,:) = handles.hhmm{rk}.evidence;
        depths{rk} = 0;
        k{rk} = 1;
    end
end
evidenceprompt(evidence,depths,k);


% --- Executes on selection change in hist2dprompts.
function hist2dprompts_Callback(hObject, eventdata, handles)
% hObject    handle to hist2dprompts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns hist2dprompts contents as cell array
%        contents{get(hObject,'Value')} returns selected item from hist2dprompts

stri = get(hObject,'String');
val = get(hObject,'Value');

n = handles.hhmm_current;
N = handles.tracestats.tracenum;
k = handles.hhmm{n}.inputoptions.k;
kp = handles.hhmm{n}.inputoptions.production;
if ~isfield(handles.hhmm{n},'pathstate')
    for rk = 1:N
        for rrk = 0:k-1
            if handles.hhmm{n}.conprob(rk,rrk+1)==max(handles.hhmm{n}.conprob(rk,:))
                for zk = 1:kp
                    handles.hhmm{n}.pathstate{rk}(handles.hhmm{n}.vbem{rk}.ideals==handles.hhmm{n}.vbem{rk}.m(zk)) = zk+rrk*kp;
                    handles.hhmm{n}.pathstate{rk} = handles.hhmm{n}.pathstate{rk}(:);
                end
            end
        end
    end
end
% compatibility
paths = cell(N,1);
switch stri{val};
    case 'Upper graph'
        for rk = 1:N
            if size(handles.hhmm{n}.pathstate{rk},2)>1
                paths{rk} = handles.hhmm{n}.pathstate{rk}(:,end);
            else
                paths{rk} = handles.hhmm{n}.pathstate{rk};
            end
        end
        dwells = hhmm_getdwells(paths);
        
        count = 0;
        for rk = 1:size(dwells,1)
            if rk ~= 1
                if dwells(rk,1)~=dwells(rk-1,1)
                    count = 0;
                end
            end
            loc(rk,:) = [dwells(rk,1) count];
            count = count+dwells(rk,4);
        end
        out = hist2Dprompt;
        start = out.bottom;
        finish = out.top;
        a = loc((dwells(:,2)==start),:);
        b = loc((dwells(:,3)==finish),:);
        limits = [handles.plotstuff.yaxis 10];
        axis(handles.graph1);
        axis on
        try
            [map,t] = fxn_hist2D(handles.hhmm{n}.data,paths,a,b,limits);
        catch
            cla
            t = 1:20;
            map = zeros(20,20);
        end
        kernel = fspecial('gaussian',5,2);
        smap = imfilter(map,kernel,'same');
        pcolor(handles.graph1,((1:size(smap,2))-11),t,smap./max(smap(:,10)));
        shading(handles.graph1,'interp');
        caxis([.01/2,1]);
        cmap = colormap('jet');
        for rk = 1:1
            cmap(rk,:) = [252 250 196]/256;
        end
        for rk = 2:2
            cmap(rk,:) = [1 1 1];
        end
        colormap(cmap);
        h = colorbar(handles.graph1);
        set(h,'ylim',[0 1]);
        set(handles.graph1,'xlim',[-10 handles.plotstuff.hxaxis(end)],'xscale','linear')
    case 'Lower graph'
        for rk = 1:N
            if size(handles.hhmm{n}.pathstate{rk},2)>1
                paths{rk} = handles.hhmm{n}.pathstate{rk}(:,end);
            else
                paths{rk} = handles.hhmm{n}.pathstate{rk};
            end
        end
        dwells = hhmm_getdwells(paths);
        
        count = 0;
        for rk = 1:size(dwells,1)
            if rk ~= 1
                if dwells(rk,1)~=dwells(rk-1,1)
                    count = 0;
                end
            end
            loc(rk,:) = [dwells(rk,1) count];
            count = count+dwells(rk,4);
        end
        out = hist2Dprompt;
        start = out.bottom;
        finish = out.top;
        a = loc((dwells(:,2)==start),:);
        b = loc((dwells(:,3)==finish),:);
        limits = [handles.plotstuff.yaxis 10];
        axis(handles.graph1);
        axis on
        cla
        try
            [map,t] = fxn_hist2D(handles.hhmm{n}.data,paths,a,b,limits);
        catch
            t = 1:20;
            map = zeros(20,20);
        end
        kernel = fspecial('gaussian',5,2);
        smap = imfilter(map,kernel,'same');
        pcolor(handles.graph2,((1:size(smap,2))-11),t,smap./max(smap(:,10)));
        shading(handles.graph2,'interp');
        caxis([.01/2,1]);
        cmap = colormap('jet');
        for rk = 1:1
            cmap(rk,:) = [252 250 196]/256;
        end
        for rk = 2:2
            cmap(rk,:) = [1 1 1];
        end
        colormap(cmap)
        hh = colorbar(handles.graph2);
        set(hh,'ylim',[0 1]);
        set(handles.graph2,'xlim',[-10 handles.plotstuff.hxaxis(end)],'xscale','linear')
end
axis(handles.modelscheme);
axis off;

% --- Executes during object creation, after setting all properties.
function hist2dprompts_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hist2dprompts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in dimensions.
function dimensions_Callback(hObject, eventdata, handles)
% hObject    handle to dimensions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    handles.plotstuff.coordlog = ~handles.plotstuff.coordlog;
catch
    handles.plotstuff.coordlog = 1==1;
end
b = handles.hhmm_current;
handles.plotstuff.instructions = 'coordplot';

try
    handles.modelscheme;
    cla;
    plotdots(handles.modelscheme,[handles.hhmm{b}.inputoptions.depth_vector handles.hhmm{b}.inputoptions.production],handles.hhmm{b}.inputoptions.k);
% make sure everything's the same if possible
catch
    handles.modelscheme;
    cla;
end

if handles.plotstuff.coordlog
    out = coordinateprompt([handles.hhmm{b}.inputoptions.depth_vector handles.hhmm{b}.inputoptions.production],handles.plotstuff.instructions);
    handles.plotstuff.coordlevel = out;
    if strcmp(out,'Direct')
        temp = [handles.hhmm{b}.inputoptions.depth_vector handles.hhmm{b}.inputoptions.production];
        handles.plotstuff.coordcolors = gray(temp(end)*handles.hhmm{b}.inputoptions.k+4);
        handles.plotstuff.coordcolors = handles.plotstuff.coordcolors(3:end-2,:);
        colordots(handles.modelscheme,handles.plotstuff.coordcolors,temp,handles.hhmm{b}.inputoptions.k,length(temp),1==1)
    else
        temp = [handles.hhmm{b}.inputoptions.depth_vector handles.hhmm{b}.inputoptions.production];
        handles.plotstuff.coordcolors = gray(prod(temp(1:out))*handles.hhmm{b}.inputoptions.k+4);
        handles.plotstuff.coordcolors = handles.plotstuff.coordcolors(3:end-2,:);
        colordots(handles.modelscheme,handles.plotstuff.coordcolors,temp,handles.hhmm{b}.inputoptions.k,out,1==0)
    end
end
guidata(hObject,handles);
scopeplot(handles);


% --- Executes on button press in xaxisprompt.
function xaxisprompt_Callback(hObject, eventdata, handles)
% hObject    handle to xaxisprompt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

out = yaxisprompt(1);
handles.plotstuff.xaxis = [1 out.top];
total = length(handles.tracestats.data.colors(1).traces(:,1));
steps = out.top./(total)/20;
set(handles.xaxisslider,'sliderstep',[steps steps],'max',1,'min',out.top./total,'Value',out.top./total);
set(handles.timetext,'string','Pan time axis');
guidata(hObject,handles);
scopeplot(handles);


% --- Executes on slider movement.
function xaxisslider_Callback(hObject, eventdata, handles)
% hObject    handle to xaxisslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

a = get(hObject,'Value');
total = length(handles.tracestats.data.colors(1).traces(:,1));
handles.plotstuff.xaxis = [a*total-diff(handles.plotstuff.xaxis) a*total];
guidata(hObject,handles);
scopeplot(handles);

% --- Executes during object creation, after setting all properties.
function xaxisslider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xaxisslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in hist2dxaxis.
function hist2dxaxis_Callback(hObject, eventdata, handles)
% hObject    handle to hist2dxaxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

out = yaxisprompt(1);
handles.plotstuff.hxaxis = out.top;
set(handles.graph1,'xlim',[-10 handles.plotstuff.hxaxis(end)],'xscale','linear')
set(handles.graph2,'xlim',[-10 handles.plotstuff.hxaxis(end)],'xscale','linear')

guidata(hObject,handles);


% --- Executes on button press in loadbutton.
function loadbutton_Callback(hObject, eventdata, handles)
% hObject    handle to loadbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fn,pn,~] = uigetfile('*.mat');
n = handles.hhmm_current;
in = load([pn fn]);
handles.hhmm{n} = in.out;
disp('Loaded!');
plotdots(handles.modelscheme,handles.hhmm{n}.depth_vector,handles.hhmm{n}.inputoptions.k);
guidata(hObject,handles);


% --- Executes on button press in parambutton.
function parambutton_Callback(hObject, eventdata, handles)
% hObject    handle to parambutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
