function varargout = settings(varargin)
% SETTINGS MATLAB code for settings.fig
%      SETTINGS, by itself, creates a new SETTINGS or raises the existing
%      singleton*.
%
%      H = SETTINGS returns the handle to a new SETTINGS or the handle to
%      the existing singleton*.
%
%      SETTINGS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SETTINGS.M with the given input arguments.
%
%      SETTINGS('Property','Value',...) creates a new SETTINGS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before settings_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to settings_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help settings

% Last Modified by GUIDE v2.5 24-Dec-2017 22:45:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @settings_OpeningFcn, ...
                   'gui_OutputFcn',  @settings_OutputFcn, ...
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


% --- Executes just before settings is made visible.
function settings_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to settings (see VARARGIN)

% Choose default command line output for settings
handles.output = hObject;

 handles.filename = '';
 handles.pathname = '';
 handles.padding ='';
 handles.lambda = '';
 handles.output_sigma_factor = '';
 handles.defaultROI = '';
 handles.isIP = '';
 handles.opticalFlow = '';
 handles.kernel = '';
 handles.features = '';
 handles.cell_size = '';
 handles.sent_info = 0;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes settings wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = settings_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.padding;   
varargout{2} = handles.lambda;
varargout{3} = handles.output_sigma_factor;
varargout{4} = handles.defaultROI;
varargout{5} = handles.isIP;
varargout{6} = handles.opticalFlow;
varargout{7} = handles.kernel;
varargout{8} = handles.pathname;
varargout{9} = handles.filename;
varargout{10} = handles.features;
varargout{11} = handles.cell_size;
delete(handles.figure1);


function edit_padding_Callback(hObject, eventdata, handles)
% hObject    handle to edit_padding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_padding as text
%        str2double(get(hObject,'String')) returns contents of edit_padding as a double


% --- Executes during object creation, after setting all properties.
function edit_padding_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_padding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_lambda_Callback(hObject, eventdata, handles)
% hObject    handle to edit_lambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_lambda as text
%        str2double(get(hObject,'String')) returns contents of edit_lambda as a double


% --- Executes during object creation, after setting all properties.
function edit_lambda_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_lambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_output_sigma_factor_Callback(hObject, eventdata, handles)
% hObject    handle to edit_output_sigma_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_output_sigma_factor as text
%        str2double(get(hObject,'String')) returns contents of edit_output_sigma_factor as a double


% --- Executes during object creation, after setting all properties.
function edit_output_sigma_factor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_output_sigma_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_ROI.
function checkbox_ROI_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_ROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_ROI


% --- Executes on button press in checkbox_isIP.
function checkbox_isIP_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_isIP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_isIP


% --- Executes on selection change in popupmenu_of.
function popupmenu_of_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_of (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_of contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_of
of_type = get(handles.popupmenu_of,'Value');
posHS = get(handles.uipanel_HS,'Position');
posLK = get(handles.uipanel_LK,'Position');
if of_type == 2
    set(handles.uipanel_HS,'visible','off');
    set(handles.uipanel_HS,'Position',posLK);
    set(handles.uipanel_LK,'Position',posHS);
    set(handles.uipanel_LK,'visible','on');
else
    set(handles.uipanel_LK,'visible','off');
    set(handles.uipanel_HS,'Position',posLK);
    set(handles.uipanel_LK,'Position',posHS);
    set(handles.uipanel_HS,'visible','on');
end;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_of_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_of (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_maxIteration_Callback(hObject, eventdata, handles)
% hObject    handle to edit_maxIteration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_maxIteration as text
%        str2double(get(hObject,'String')) returns contents of edit_maxIteration as a double


% --- Executes during object creation, after setting all properties.
function edit_maxIteration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_maxIteration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_smoothness_Callback(hObject, eventdata, handles)
% hObject    handle to edit_smoothness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_smoothness as text
%        str2double(get(hObject,'String')) returns contents of edit_smoothness as a double


% --- Executes during object creation, after setting all properties.
function edit_smoothness_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_smoothness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_noiseThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to edit_noiseThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_noiseThreshold as text
%        str2double(get(hObject,'String')) returns contents of edit_noiseThreshold as a double


% --- Executes during object creation, after setting all properties.
function edit_noiseThreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_noiseThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_kernel.
function popupmenu_kernel_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_kernel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_kernel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_kernel
pos_left = [3.2, 1.2307692307692308, 65.6, 14.0];
pos_right = [73.2, 1.2307692307692308, 65.6, 14.0];

kernel_type = get(handles.popupmenu_kernel,'Value');
feature_type = get(handles.popupmenu_feature,'Value');
if kernel_type ==1 
    set(handles.uipanel_kernel,'Visible','on');
    
    set(handles.uipanel_gauss,'Visible','on');
    set(handles.uipanel_gauss,'Position',pos_left);
    
    set(handles.uipanel_poly,'Visible','off');
    set(handles.uipanel_poly,'Position',pos_right);
    
    if feature_type == 1
        set(handles.edit_gauss_interp,'String','0.02');
        set(handles.edit_sigma,'String','0.5');
        set(handles.edit_gauss_cell,'String','4');
    else
        set(handles.edit_gauss_interp,'String','0.075');
        set(handles.edit_sigma,'String','0.2');
        set(handles.edit_gauss_cell,'String','1');
    end              
elseif kernel_type == 2
    set(handles.uipanel_kernel,'Visible','on');
    
    set(handles.uipanel_poly,'Visible','on');
    set(handles.uipanel_poly,'Position',pos_left);
    
    set(handles.uipanel_gauss,'Visible','off');
    set(handles.uipanel_gauss,'Position',pos_right);
    
    if feature_type == 1
        set(handles.edit_poly_interp,'String','0.02');
        set(handles.edit_poly_a,'String','1');
        set(handles.edit_poly_b,'String','7');
        set(handles.edit_poly_cell,'String','4');
    else
        set(handles.edit_poly_interp,'String','0.075');
        set(handles.edit_poly_a,'String','1');
        set(handles.edit_poly_b,'String','9');
        set(handles.edit_poly_cell,'String','1');
    end   
else
    set(handles.uipanel_kernel,'Visible','off');
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_kernel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_kernel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_feature.
function popupmenu_feature_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_feature (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_feature contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_feature
feature_type = get(handles.popupmenu_feature,'Value');
if feature_type == 1
    set(handles.uipanel_feature,'visible','on');
else
    set(handles.uipanel_feature,'visible','off');
end;
pos_left = [3.2, 1.2307692307692308, 65.6, 14.0];
pos_right = [73.2, 1.2307692307692308, 65.6, 14.0];

kernel_type = get(handles.popupmenu_kernel,'Value');
if kernel_type ==1 
    set(handles.uipanel_kernel,'Visible','on');
    
    set(handles.uipanel_gauss,'Visible','on');
    set(handles.uipanel_gauss,'Position',pos_left);
    
    set(handles.uipanel_poly,'Visible','off');
    set(handles.uipanel_poly,'Position',pos_right);
    
    if feature_type == 1
        set(handles.edit_gauss_interp,'String','0.02');
        set(handles.edit_sigma,'String','0.5');
        set(handles.edit_gauss_cell,'String','4');
    else
        set(handles.edit_gauss_interp,'String','0.075');
        set(handles.edit_sigma,'String','0.2');
        set(handles.edit_gauss_cell,'String','1');
    end              
elseif kernel_type == 2
    set(handles.uipanel_kernel,'Visible','on');
    
    set(handles.uipanel_poly,'Visible','on');
    set(handles.uipanel_poly,'Position',pos_left);
    
    set(handles.uipanel_gauss,'Visible','off');
    set(handles.uipanel_gauss,'Position',pos_right);
    
    if feature_type == 1
        set(handles.edit_poly_interp,'String','0.02');
        set(handles.edit_poly_a,'String','1');
        set(handles.edit_poly_b,'String','7');
        set(handles.edit_poly_cell,'String','4');
    else
        set(handles.edit_poly_interp,'String','0.075');
        set(handles.edit_poly_a,'String','1');
        set(handles.edit_poly_b,'String','9');
        set(handles.edit_poly_cell,'String','1');
    end   
else
    set(handles.uipanel_kernel,'Visible','off');
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_feature_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_feature (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_gauss_interp_Callback(hObject, eventdata, handles)
% hObject    handle to edit_gauss_interp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_gauss_interp as text
%        str2double(get(hObject,'String')) returns contents of edit_gauss_interp as a double


% --- Executes during object creation, after setting all properties.
function edit_gauss_interp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_gauss_interp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_sigma_Callback(hObject, eventdata, handles)
% hObject    handle to edit_sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_sigma as text
%        str2double(get(hObject,'String')) returns contents of edit_sigma as a double


% --- Executes during object creation, after setting all properties.
function edit_sigma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double


% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_hog_Callback(hObject, eventdata, handles)
% hObject    handle to edit_hog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_hog as text
%        str2double(get(hObject,'String')) returns contents of edit_hog as a double


% --- Executes during object creation, after setting all properties.
function edit_hog_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_hog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit19_Callback(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit19 as text
%        str2double(get(hObject,'String')) returns contents of edit19 as a double


% --- Executes during object creation, after setting all properties.
function edit19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit20_Callback(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit20 as text
%        str2double(get(hObject,'String')) returns contents of edit20 as a double


% --- Executes during object creation, after setting all properties.
function edit20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_gauss_cell_Callback(hObject, eventdata, handles)
% hObject    handle to edit_gauss_cell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_gauss_cell as text
%        str2double(get(hObject,'String')) returns contents of edit_gauss_cell as a double


% --- Executes during object creation, after setting all properties.
function edit_gauss_cell_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_gauss_cell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_poly_a_Callback(hObject, eventdata, handles)
% hObject    handle to edit_poly_a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_poly_a as text
%        str2double(get(hObject,'String')) returns contents of edit_poly_a as a double


% --- Executes during object creation, after setting all properties.
function edit_poly_a_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_poly_a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_poly_b_Callback(hObject, eventdata, handles)
% hObject    handle to edit_poly_b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_poly_b as text
%        str2double(get(hObject,'String')) returns contents of edit_poly_b as a double


% --- Executes during object creation, after setting all properties.
function edit_poly_b_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_poly_b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_poly_interp_Callback(hObject, eventdata, handles)
% hObject    handle to edit_poly_interp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_poly_interp as text
%        str2double(get(hObject,'String')) returns contents of edit_poly_interp as a double


% --- Executes during object creation, after setting all properties.
function edit_poly_interp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_poly_interp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_poly_cell_Callback(hObject, eventdata, handles)
% hObject    handle to edit_poly_cell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_poly_cell as text
%        str2double(get(hObject,'String')) returns contents of edit_poly_cell as a double


% --- Executes during object creation, after setting all properties.
function edit_poly_cell_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_poly_cell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_load.
function pushbutton_load_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( ...
   {'*.avi;*.mpg;*.mpeg;*.mp4;*.mkv','Video Files (*.avi,*.mpg,*.mpeg,*.mp4,*.mkv)';
     '*.*',  'All Files (*.*)'}, ...
     'Select a video file');
 fullname = strcat(pathname, filename);
 set(handles.text_file,'String',fullname);
 set(handles.pushbutton_last,'Visible','on');
 handles.filename = filename;
 handles.pathname = pathname;
guidata(hObject,handles);

% --- Executes on button press in pushbutton_last.
function pushbutton_last_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_last (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
padding = str2double(get(handles.edit_padding,'String'));
lambda = str2double(get(handles.edit_lambda,'String'));
output_sigma_factor = str2double(get(handles.edit_output_sigma_factor,'String'));
defaultROI = get(handles.checkbox_ROI,'Value');
isIP = get(handles.checkbox_isIP,'Value');

oftype = get(handles.popupmenu_of,'Value');
if oftype == 1
    opticalFlow.type = 'Horn-Schunck';
    opticalFlow.maxIteration = str2double(get(handles.edit_maxIteration,'String'));
    opticalFlow.smoothness = str2double(get(handles.edit_smoothness,'String'));
else
    opticalFlow.type = 'Lucas-Kanade';
    opticalFlow.noiseThreshold = str2double(get(handles.edit_noiseThreshold,'String'));
end

kerneltype = get(handles.popupmenu_kernel,'Value');
if kerneltype == 1
    kernel.type = 'gaussian';
    kernel.sigma = str2double(get(handles.edit_sigma,'String'));
    kernel.interp_factor = str2double(get(handles.edit_gauss_interp,'String'));
    cell_size = str2double(get(handles.edit_gauss_cell,'String'));
elseif kerneltype == 2
    kernel.type = 'gaussian';
    kernel.poly_a = str2double(get(handles.edit_poly_a,'String'));
	kernel.poly_b = str2double(get(handles.edit_poly_b,'String'));
    kernel.interp_factor = str2double(get(handles.edit_poly_interp,'String'));
    cell_size = str2double(get(handles.edit_poly_cell,'String'));
else
    kernel.type = 'linear';
end

featuretype = get(handles.popupmenu_feature,'Value');
if featuretype == 1
    features.gray = false;
    features.hog = true;
    features.hog_orientations = str2double(get(handles.edit_hog,'String'));
else
    features.gray = true;
    features.hog = false;
end

%% check
errorFlag = false;
errorString = '';

if isnan(padding)
    errorFlag = true;
    errorString = 'Padding ';
end

if isnan(lambda)
    errorFlag = true;
    errorString = strcat(errorString,' Lambda ');
end

if isnan(output_sigma_factor)
    errorFlag = true;
    errorString = strcat(errorString,' Output Sigma Factor ');
end

if oftype == 1
    if isnan(opticalFlow.maxIteration)
        errorFlag = true;
        errorString = strcat(errorString,' Max Iteration ');
    end
    if isnan(opticalFlow.smoothness)
        errorFlag = true;
        errorString = strcat(errorString,' Smoothness ');
    end
else
    if isnan(opticalFlow.noiseThreshold)
        errorFlag = true;
        errorString = strcat(errorString,' Noise Threshold ');
    end
end

if featuretype == 1
    if isnan(features.hog_orientations)
        errorFlag = true;
        errorString = strcat(errorString,' Hog Orientations ');
    end
end

if isnan(kernel.interp_factor)
    errorFlag = true;
    errorString = strcat(errorString,' Sigma ');
end
   
if isnan(cell_size)
    errorFlag = true;
    errorString = strcat(errorString,' Cell size ');
end

if kerneltype == 1
    if isnan(kernel.sigma)
        errorFlag = true;
        errorString = strcat(errorString,' Sigma ');
    end
elseif kerneltype == 2
    if isnan(kernel.poly_a)
        errorFlag = true;
        errorString = strcat(errorString,' Additive Term ');
    end
    if isnan(kernel.poly_b)
        errorFlag = true;
        errorString = strcat(errorString,' Exponent ');
    end    
end

if strcmp(handles.filename,'') || strcmp(handles.pathname,'')
    errorFlag = true;
    errorString = strcat(errorString,' File Selection ');
end

if errorFlag
    errorString = strcat(errorString,' is not selected, entered or have incompatible data. Please check!');
    msgbox(errorString, 'Error','error');
else
    handles.padding = padding;
    handles.lambda = lambda;
    handles.output_sigma_factor = output_sigma_factor;
    handles.defaultROI = defaultROI;
    handles.isIP = isIP;
    handles.opticalFlow = opticalFlow;
    handles.kernel = kernel;
    handles.features = features;
    handles.cell_size = cell_size;
     handles.sent_info = 1;

    msgbox('Data is written, you can exit!', 'Successful','warn');
end
guidata(hObject,handles);



% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 if ~handles.sent_info
    button = questdlg('You have not sent the settings information. This can cause some errors. Are you sre to exit?','Warning','Exit','Return to settings','Return to settings');
    if strcmp(button,'Return to settings')
        return;
    end
 end
 
% Hint: delete(hObject) closes the figure
if isequal(get(hObject, 'waitstatus'), 'waiting')
    uiresume(hObject);
else
    delete(hObject);
end
