function varargout = guiFracTend(varargin)
% GUIFRACTEND MATLAB code for guiFracTend.fig
%      GUIFRACTEND, by itself, creates a new GUIFRACTEND or raises the existing
%      singleton*.
%
%      H = GUIFRACTEND returns the handle to a new GUIFRACTEND or the handle to
%      the existing singleton*.
%
%      GUIFRACTEND('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUIFRACTEND.M with the given input arguments.
%
%      GUIFRACTEND('Property','Value',...) creates a new GUIFRACTEND or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before guiFracTend_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to guiFracTend_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help guiFracTend

% Last Modified by GUIDE v2.5 10-Jul-2018 09:40:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @guiFracTend_OpeningFcn, ...
                   'gui_OutputFcn',  @guiFracTend_OutputFcn, ...
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


% --- Executes just before guiFracTend is made visible.
function guiFracTend_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to guiFracTend (see VARARGIN)

% Choose default command line output for guiFracTend
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes guiFracTend wait for user response (see UIRESUME)
% uiwait(handles.figure1);
handles.FracTendversion = '1.0' ; 
disp(['FracTend version ', handles.FracTendversion]) ; 
% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = guiFracTend_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pbRun.
function pbRun_Callback(hObject, eventdata, handles)
% hObject    handle to pbRun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

flagError = false ; 

sValue = get(handles.editMu, 'String') ; 
if isnan(str2double(sValue)) || str2double(sValue) < 0 || str2double(sValue) > 2 
    hError = errordlg('Friction coefficient must be a positive decimal (e.g. 0.1-1.0)', ...
                        'Input error', 'modal') ; 
    uicontrol(handles.editMu) ; 
    flagError = true ; 
    return ; 
end

sValue = get(handles.editC0, 'String') ; 
if isnan(str2double(sValue)) || str2double(sValue) < 0  
    hError = errordlg('Cohesion must be a positive number (e.g. 0.0-200 MPa)', ...
                        'Input error', 'modal') ; 
    uicontrol(handles.editC0) ; 
    flagError = true ; 
    return ; 
end

sValue = get(handles.editS1, 'String') ; 
if isnan(str2double(sValue))  
    hError = errordlg('Sigma 1 (most compressive) must be a decimal number (e.g. 0.0-200 MPa)', ...
                        'Input error', 'modal') ; 
    uicontrol(handles.editS1) ; 
    flagError = true ; 
    return ; 
end

sValue = get(handles.editS2, 'String') ; 
if isnan(str2double(sValue))  
    hError = errordlg('Sigma 2 (intermediate stress) must be a decimal number (e.g. 0.0-200 MPa)', ...
                        'Input error', 'modal') ; 
    uicontrol(handles.editS2) ; 
    flagError = true ; 
    return ; 
end

sValue = get(handles.editS3, 'String') ; 
if isnan(str2double(sValue)) 
    hError = errordlg('Sigma 3 (least compressive) must be a decimal number (e.g. 0.0-200 MPa)', ...
                        'Input error', 'modal') ; 
    uicontrol(handles.editS3) ; 
    flagError = true ; 
    return ; 
end

sValue = get(handles.editPf, 'String') ; 
if isnan(str2double(sValue)) 
    hError = errordlg('Pf (pore fluid pressure) must be a decimal number (e.g. 0.0-50 MPa)', ...
                        'Input error', 'modal') ; 
    uicontrol(handles.editPf) ; 
    flagError = true ; 
    return ; 
end

sValueS1 = get(handles.editS1, 'String') ; 
sValueS2 = get(handles.editS2, 'String') ; 
sValueS3 = get(handles.editS3, 'String') ; 
if str2double(sValueS1) < str2double(sValueS2)     
    hError = errordlg('Sigma 1 (most compressive) must be greater than Sigma 2 (intermediate)', ...
                        'Input error', 'modal') ; 
    uicontrol(handles.editS1) ; 
    flagError = true ; 
    return ; 
end

if str2double(sValueS2) < str2double(sValueS3)
    hError = errordlg('Sigma 2 (intermediate) must be greater than Sigma 3 (least compressive)', ...
                        'Input error', 'modal') ; 
    uicontrol(handles.editS2) ; 
    flagError = true ; 
    return ; 
end 

if str2double(sValueS1) < str2double(sValueS3)     
    hError = errordlg('Sigma 1 (most compressive) must be greater than Sigma 3 (least compressive)', ...
                        'Input error', 'modal') ; 
    uicontrol(handles.editS1) ; 
    flagError = true ; 
    return ; 
end

sValue = get(handles.editS1Trend, 'String') ; 
if isnan(str2double(sValue)) || str2double(sValue) < 0 || str2double(sValue) > 360 
    hError = errordlg('Sigma 1 trend must be a (positive) azimuth (e.g. 0-360 degrees)', ...
                        'Input error', 'modal') ; 
    uicontrol(handles.editS1Trend) ; 
    flagError = true ; 
    return ; 
end

sValue = get(handles.editS1Plunge, 'String') ; 
if isnan(str2double(sValue)) || str2double(sValue) < 0 || str2double(sValue) > 90 
    hError = errordlg('Sigma 1 trend must be a (positive) plunge (e.g. 0-90 degrees)', ...
                        'Input error', 'modal') ; 
    uicontrol(handles.editS1Plunge) ; 
    flagError = true ; 
    return ; 
end

sValue = get(handles.editS3Trend, 'String') ; 
if isnan(str2double(sValue)) || str2double(sValue) < 0 || str2double(sValue) > 360 
    hError = errordlg('Sigma 3 trend must be a (positive) azimuth (e.g. 0-360 degrees)', ...
                        'Input error', 'modal') ; 
    uicontrol(handles.editS3Trend) ; 
    flagError = true ; 
    return ; 
end

if ~flagError

    %   call the various maps & graphs 
    flag_TsStereo = get(handles.checkboxTsStereo, 'Value') ; 
    flag_TdStereo = get(handles.checkboxTdStereo, 'Value') ; 
    flag_SfStereo = get(handles.checkboxSfStereo, 'Value') ; 
    flag_OAStereo = get(handles.checkboxOAStereo, 'Value') ; 

    flag_TsMohr = get(handles.checkboxTsMohr, 'Value') ; 
    flag_TdMohr = get(handles.checkboxTdMohr, 'Value') ; 
    flag_SfMohr = get(handles.checkboxSfMohr, 'Value') ; 
    flag_OAMohr = get(handles.checkboxOAMohr, 'Value') ; 
    
    %   get value from drop down list 
    list = get(handles.popupIncrement, 'String') ; 
    nIncrement = str2double(list{get(handles.popupIncrement, 'Value')}) ;

    %   get values from edit text boxes 
    nS1 = str2double(get(handles.editS1, 'String')) ;
    nS2 = str2double(get(handles.editS2, 'String')) ;
    nS3 = str2double(get(handles.editS3, 'String')) ;
    nPf = str2double(get(handles.editPf, 'String')) ;
    nS1Trend = str2double(get(handles.editS1Trend, 'String')) ;
    nS1Plunge = str2double(get(handles.editS1Plunge, 'String')) ;
    nS3Trend = str2double(get(handles.editS3Trend, 'String')) ;
    nMu = str2double(get(handles.editMu, 'String')) ;
    nC0 = str2double(get(handles.editC0, 'String')) ;
    
    runFracTend(handles.selfile, handles.selpath, ... 
                nS1, nS2, nS3, nPf, nS1Trend, nS1Plunge, nS3Trend, nMu, nC0, ... 
                flag_TsStereo, flag_TdStereo, flag_SfStereo, flag_OAStereo, ...
                flag_TsMohr, flag_TdMohr, flag_SfMohr, flag_OAMohr) ; 
    
end 

% --- Executes on button press in pbExit.
function pbExit_Callback(hObject, eventdata, handles)
% hObject    handle to pbExit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close all ; 
close(gcf) ; 


function editFilename_Callback(hObject, eventdata, handles)
% hObject    handle to editFilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editFilename as text
%        str2double(get(hObject,'String')) returns contents of editFilename as a double


% --- Executes during object creation, after setting all properties.
function editFilename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbBrowse.
function pbBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to pbBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[ handles.selfile, handles.selpath ] = uigetfile('*.*', 'Select input file(s) for FracTend') ; 

if ~isempty(handles.selfile)
    
    sFileName = char(handles.selfile) ; 
    set(handles.editFilename, 'String', sFileName) ; 
    set(handles.pbRun, 'Enable', 'on') ; 

else
    
    set(handles.editFilename, 'String', '(no file selected)') ; 
    set(handles.pbRun, 'Enable', 'off') ; 
    
end

% Update handles structure
guidata(hObject, handles);

% --- Executes on selection change in popupIncrement.
function popupIncrement_Callback(hObject, eventdata, handles)
% hObject    handle to popupIncrement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupIncrement contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupIncrement


% --- Executes during object creation, after setting all properties.
function popupIncrement_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupIncrement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editMu_Callback(hObject, eventdata, handles)
% hObject    handle to editMu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMu as text
%        str2double(get(hObject,'String')) returns contents of editMu as a double


% --- Executes during object creation, after setting all properties.
function editMu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editC0_Callback(hObject, eventdata, handles)
% hObject    handle to editC0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editC0 as text
%        str2double(get(hObject,'String')) returns contents of editC0 as a double


% --- Executes during object creation, after setting all properties.
function editC0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editC0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editS1_Callback(hObject, eventdata, handles)
% hObject    handle to editS1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editS1 as text
%        str2double(get(hObject,'String')) returns contents of editS1 as a double


% --- Executes during object creation, after setting all properties.
function editS1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editS1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editS1Trend_Callback(hObject, eventdata, handles)
% hObject    handle to editS1Trend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editS1Trend as text
%        str2double(get(hObject,'String')) returns contents of editS1Trend as a double


% --- Executes during object creation, after setting all properties.
function editS1Trend_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editS1Trend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editS1Plunge_Callback(hObject, eventdata, handles)
% hObject    handle to editS1Plunge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editS1Plunge as text
%        str2double(get(hObject,'String')) returns contents of editS1Plunge as a double


% --- Executes during object creation, after setting all properties.
function editS1Plunge_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editS1Plunge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editS2_Callback(hObject, eventdata, handles)
% hObject    handle to editS2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editS2 as text
%        str2double(get(hObject,'String')) returns contents of editS2 as a double


% --- Executes during object creation, after setting all properties.
function editS2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editS2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editS3_Callback(hObject, eventdata, handles)
% hObject    handle to editS3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editS3 as text
%        str2double(get(hObject,'String')) returns contents of editS3 as a double


% --- Executes during object creation, after setting all properties.
function editS3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editS3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editS3Trend_Callback(hObject, eventdata, handles)
% hObject    handle to editS3Trend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editS3Trend as text
%        str2double(get(hObject,'String')) returns contents of editS3Trend as a double


% --- Executes during object creation, after setting all properties.
function editS3Trend_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editS3Trend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editPf_Callback(hObject, eventdata, handles)
% hObject    handle to editPf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editPf as text
%        str2double(get(hObject,'String')) returns contents of editPf as a double


% --- Executes during object creation, after setting all properties.
function editPf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxTsMohr.
function checkboxTsMohr_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxTsMohr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxTsMohr


% --- Executes on button press in checkboxTdMohr.
function checkboxTdMohr_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxTdMohr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxTdMohr


% --- Executes on button press in checkboxSfMohr.
function checkboxSfMohr_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxSfMohr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxSfMohr


% --- Executes on button press in checkboxOAMohr.
function checkboxOAMohr_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxOAMohr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxOAMohr


% --- Executes on button press in checkboxTsStereo.
function checkboxTsStereo_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxTsStereo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxTsStereo


% --- Executes on button press in checkboxTdStereo.
function checkboxTdStereo_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxTdStereo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxTdStereo


% --- Executes on button press in checkboxSfStereo.
function checkboxSfStereo_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxSfStereo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxSfStereo


% --- Executes on button press in checkboxOAStereo.
function checkboxOAStereo_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxOAStereo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxOAStereo
