function varargout = evidenceprompt(varargin)
% EVIDENCEPROMPT MATLAB code for evidenceprompt.fig
%      EVIDENCEPROMPT, by itself, creates a new EVIDENCEPROMPT or raises the existing
%      singleton*.
%
%      H = EVIDENCEPROMPT returns the handle to a new EVIDENCEPROMPT or the handle to
%      the existing singleton*.
%
%      EVIDENCEPROMPT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EVIDENCEPROMPT.M with the given input arguments.
%
%      EVIDENCEPROMPT('Property','Value',...) creates a new EVIDENCEPROMPT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before evidenceprompt_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to evidenceprompt_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help evidenceprompt

% Last Modified by GUIDE v2.5 28-Mar-2018 11:34:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @evidenceprompt_OpeningFcn, ...
                   'gui_OutputFcn',  @evidenceprompt_OutputFcn, ...
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


% --- Executes just before evidenceprompt is made visible.
function evidenceprompt_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to evidenceprompt (see VARARGIN)

% Choose default command line output for evidenceprompt
handles.output = hObject;

% Update handles structure

evidence = varargin{1};
depths = varargin{2};
k = varargin{3};
for rk = 1:length(depths)
    lendep(rk,:) = length(depths{rk})+1;
end
axis(handles.evidence);
cla
s = nanstd(evidence(~isinf(evidence)));
m = min(evidence(~isinf(evidence)&~isnan(evidence)));
mm = max(evidence(~isinf(evidence)&~isnan(evidence)));
bar(handles.evidence,evidence,'k');
set(handles.evidence,'ylim',([m-s mm+s]),'xlim',[0.5 10.55]);
set(handles.evidence,'xticklabel',[]);
ylabel(handles.evidence,'Log Evidence Lower Bound (nats)');

axis(handles.m1);
axis off
try
    plotdots(handles.m1,depths{1},max(lendep),k{1});
end

axis(handles.m2);
axis off
try
    plotdots(handles.m2,depths{2},max(lendep),k{2});
end

axis(handles.m3);
axis off
try
    plotdots(handles.m3,depths{3},max(lendep),k{3});
end

axis(handles.m4);
axis off
try
    plotdots(handles.m4,depths{4},max(lendep),k{4});
end

axis(handles.m5);
axis off
try
    plotdots(handles.m5,depths{5},max(lendep),k{5});
end

axis(handles.m6);
axis off
try
    plotdots(handles.m6,depths{6},max(lendep),k{6});
end

axis(handles.m7);
axis off
try
    plotdots(handles.m7,depths{7},max(lendep),k{7});
end

axis(handles.m8);
axis off
try
    plotdots(handles.m8,depths{8},max(lendep),k{8});
end

axis(handles.m9);
axis off
try
    plotdots(handles.m9,depths{9},max(lendep),k{9});
end

axis(handles.m10);
axis off
try
    plotdots(handles.m10,depths{10},max(lendep),k{10});
end

guidata(hObject, handles);

% UIWAIT makes evidenceprompt wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = evidenceprompt_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function plotdots(hax,depth,ytop,k)
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
                plot(hax,tt+nk*max(tt),0*tt+rk,'.k','MarkerSize',10)
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
        plot(hax,t+nk*max(t),0*(t)+1,'.r','MarkerSize',10)
    end    
    hold off
    ylim(hax,[0 ytop])
    xlim(hax,[0 (prod(depth))*2+(k-1)*(prod(depth))*2])
    axis off
    
