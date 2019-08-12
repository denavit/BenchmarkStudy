function varargout = runBenchmarkStudyGUI(varargin)
% RUNBENCHMARKSTUDYGUI MATLAB code for runBenchmarkStudyGUI.fig
%      RUNBENCHMARKSTUDYGUI, by itself, creates a new RUNBENCHMARKSTUDYGUI or raises the existing
%      singleton*.
%
%      H = RUNBENCHMARKSTUDYGUI returns the handle to a new RUNBENCHMARKSTUDYGUI or the handle to
%      the existing singleton*.
%
%      RUNBENCHMARKSTUDYGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RUNBENCHMARKSTUDYGUI.M with the given input arguments.
%
%      RUNBENCHMARKSTUDYGUI('Property','Value',...) creates a new RUNBENCHMARKSTUDYGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before runBenchmarkStudyGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to runBenchmarkStudyGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help runBenchmarkStudyGUI

% Last Modified by GUIDE v2.5 16-Mar-2013 17:35:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @runBenchmarkStudyGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @runBenchmarkStudyGUI_OutputFcn, ...
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


% --- Executes just before runBenchmarkStudyGUI is made visible.
function runBenchmarkStudyGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to runBenchmarkStudyGUI (see VARARGIN)

% Choose default command line output for runBenchmarkStudyGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes runBenchmarkStudyGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = runBenchmarkStudyGUI_OutputFcn(hObject, eventdata, handles) 
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
set(hObject,'String',study_names());


% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2
handles = guihandles(hObject);
contents = cellstr(get(hObject,'String'));
action = contents{get(hObject,'Value')};
listbox3options = options(action);
if get(handles.listbox3,'Value') > length(listbox3options)
    set(handles.listbox3,'Value',1)
end
set(handles.listbox3,'String',listbox3options);

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
set(hObject,'String',actions());

% --- Executes on selection change in listbox3.
function listbox3_Callback(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox3


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
set(hObject,'String',{''});

% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1



% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guihandles(hObject);
listbox1contents = cellstr(get(handles.listbox1,'String'));
studyName = listbox1contents{get(handles.listbox1,'Value')};
listbox2contents = cellstr(get(handles.listbox2,'String'));
action = listbox2contents{get(handles.listbox2,'Value')};
listbox3contents = cellstr(get(handles.listbox3,'String'));
option = listbox3contents{get(handles.listbox3,'Value')};
deletePrevious = get(handles.radiobutton1,'Value');
saveFigures = get(handles.radiobutton2,'Value');
runBenchmarkStudy(studyName,action,option,deletePrevious,saveFigures)


function x = study_names()
x = {'CCFT','RCFT','SRCs','SRCw','RoundHSS','RectHSS','WFs','WFw','Maleck_SA','Maleck_WA','RC'};


function x = actions()
x = {'List Results Files',...
    '--------------------------',...
    'Analysis Interaction',...
    'Figure - Analysis Interaction',...
    'Analysis Interaction (Section)',...
    'Figure - Analysis Interaction (Section)',...    
    '--------------------------',...
    'Design Interaction',...
    'Figure - Design Interaction',...
    '--------------------------',...
    'ACI Interaction',...
    'Figure - ACI Interaction',...
    'ACI Interaction (Section)',...
    'Figure - ACI Interaction (Section)'};


function x = options(action)
switch action
    case {'Design Interaction','Figure - Design Interaction'}
        x = {...
            'AISC 2016 (DA)',...
            'AISC 2016 (DMMI)',...
            'AISC 2016 (EL)',...
            'ACI 2011',...
            'ACI 2011 (No Moment Ratio Limit)',...
            'Maleck (DA)',...
            'Maleck (EL)',...
            'Scratch'};
    case 'ACI Interaction'
        x = {'a','b','c','a - no limit','b - no limit','c - no limit'};
    case 'Figure - ACI Interaction'
        x = {'a','b','c','a - no limit','b - no limit','c - no limit','a & a - no limit','b & b - no limit','c & c - no limit'};
    case {'ACI Interaction (Section)','Figure - ACI Interaction (Section)'}
        x = {'nominal','reduced'};
    otherwise
        x = {''};
end
