function varargout = parameter_setting(varargin)
% PARAMETER_SETTING MATLAB code for parameter_setting.fig
%      PARAMETER_SETTING, by itself, creates a new PARAMETER_SETTING or raises the existing
%      singleton*.
%
%      H = PARAMETER_SETTING returns the handle to a new PARAMETER_SETTING or the handle to
%      the existing singleton*.
%
%      PARAMETER_SETTING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PARAMETER_SETTING.M with the given input arguments.
%
%      PARAMETER_SETTING('Property','Value',...) creates a new PARAMETER_SETTING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before parameter_setting_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to parameter_setting_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help parameter_setting

% Last Modified by GUIDE v2.5 14-Oct-2016 11:36:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @parameter_setting_OpeningFcn, ...
                   'gui_OutputFcn',  @parameter_setting_OutputFcn, ...
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
end


% --- Executes just before parameter_setting is made visible.
function parameter_setting_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to parameter_setting (see VARARGIN)

% Choose default command line output for parameter_setting
% Choose default command line output for eyesmFISH_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Update handles structure
load(['/functions/' 'root.mat']);
load([datafolder 'data.mat']);
handles.datafolder=datafolder;

% load p parameters from data.mat in datafolder
set(handles.checkbox1,'Value',p.savetiff);
set(handles.checkbox2,'Value',p.savefinalresult);
set(handles.checkbox9,'Value',p.rawImage16);
set(handles.checkbox8,'Value',p.closeallfigure);
set(handles.popupmenu3,'Value',p.imageformat);
set(handles.edit6,'String',num2str(p.DAPI_chal));
set(handles.edit7,'String',num2str(p.FISH_chal));
if strcmp(p.gau_type,'gaussian')==1
    set(handles.popupmenu1,'Value',1);
else
    set(handles.popupmenu1,'Value',1);
end
set(handles.edit3,'String',num2str(p.gau_hsize));
set(handles.edit4,'String',num2str(p.gau_sigma));
set(handles.edit8,'String',num2str(p.init_stack));
set(handles.edit9,'String',num2str(p.final_stack));
set(handles.edit2,'String',num2str(p.back_th_num));
set(handles.edit5,'String',num2str(p.back_back_erode));
set(handles.edit10,'String',num2str(p.back_back_open));
set(handles.edit11,'String',num2str(p.back_nuc_open));
set(handles.edit12,'String',num2str(p.back_nuc_peakthresh));
set(handles.edit13,'String',num2str(p.candi_start));
set(handles.edit14,'String',num2str(p.candi_end));
set(handles.edit15,'String',num2str(p.nearby_radius_max));
set(handles.edit16,'String',num2str(p.remove_bou_back));
set(handles.edit17,'String',num2str(p.dist_radius_ratio));

guidata(hObject, handles);

% UIWAIT makes parameter_setting wait for user response (see UIRESUME)
% uiwait(handles.figure1);
end

% --- Outputs from this function are returned to the command line.
function varargout = parameter_setting_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%varargout{1} = handles.output;

% Update handles structure

% Update handles structure
guidata(hObject, handles);

% Update handles structure
load(['/functions/' 'root.mat']);
load([datafolder 'data.mat']);
%{
load(['/functions/' 'root.mat']);

p.savetiff=get(handles.checkbox1,'Value');
p.savefinalresult=get(handles.checkbox2,'Value');
p.rawImage16=get(handles.checkbox9,'Value');
p.closeallfigure=get(handles.checkbox8,'Value');
p.imageformat=get(handles.popupmenu3,'Value');
p.DAPI_chal=str2double(get(handles.edit6,'String'));
p.FISH_chal=str2double(get(handles.edit7,'String'));
if get(handles.popupmenu1,'Value')==1;
    p.gau_type = 'gaussian';
else
    p.gau_type = 'gaussian';
end
p.gau_hsize=str2double(get(handles.edit3,'String'));
p.gau_sigma=str2double(get(handles.edit4,'String'));
p.init_stack=str2double(get(handles.edit8,'String'));
p.final_stack=str2double(get(handles.edit9,'String'));
p.back_th_num=str2double(get(handles.edit2,'String'));
p.back_back_erode=str2double(get(handles.edit5,'String'));
p.back_back_open=str2double(get(handles.edit10,'String'));
p.back_nuc_open=str2double(get(handles.edit11,'String'));
p.back_nuc_peakthresh=str2double(get(handles.edit12,'String'));
p.candi_start=str2double(get(handles.edit13,'String'));
p.candi_end=str2double(get(handles.edit14,'String'));
p.nearby_radius_max=str2double(get(handles.edit15,'String'));
p.remove_bou_back=str2double(get(handles.edit16,'String'));
p.dist_radius_ratio=str2double(get(handles.edit17,'String'));

save([datafolder 'data.mat'],'p','-append');

guidata(hObject, handles);
%}
varargout{1} = p;

end

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Update handles structure
guidata(hObject, handles);
load(['/functions/' 'root.mat']);
load([datafolder 'data.mat']);

p.savetiff=get(handles.checkbox1,'Value');
p.savefinalresult=get(handles.checkbox2,'Value');
p.rawImage16=get(handles.checkbox9,'Value');
p.closeallfigure=get(handles.checkbox8,'Value');
p.imageformat=get(handles.popupmenu3,'Value');
p.DAPI_chal=str2double(get(handles.edit6,'String'));
p.FISH_chal=str2double(get(handles.edit7,'String'));
if get(handles.popupmenu1,'Value')==1;
    p.gau_type = 'gaussian';
else
    p.gau_type = 'gaussian';
end
p.gau_hsize=str2double(get(handles.edit3,'String'));
p.gau_sigma=str2double(get(handles.edit4,'String'));
p.init_stack=str2double(get(handles.edit8,'String'));
p.final_stack=str2double(get(handles.edit9,'String'));
p.back_th_num=str2double(get(handles.edit2,'String'));
p.back_back_erode=str2double(get(handles.edit5,'String'));
p.back_back_open=str2double(get(handles.edit10,'String'));
p.back_nuc_open=str2double(get(handles.edit11,'String'));
p.back_nuc_peakthresh=str2double(get(handles.edit12,'String'));
p.candi_start=str2double(get(handles.edit13,'String'));
p.candi_end=str2double(get(handles.edit14,'String'));
p.nearby_radius_max=str2double(get(handles.edit15,'String'));
p.remove_bou_back=str2double(get(handles.edit16,'String'));
p.dist_radius_ratio=str2double(get(handles.edit17,'String'));

save([datafolder 'data.mat'],'p','-append');
guidata(hObject, handles);
varargout{1} = p;
% Hint: delete(hObject) closes the figure
delete(hObject);
end

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Choose default command line output for parameter_setting
handles.output = hObject;

% Update handles structure
load(['/functions/' 'root.mat']);
load([datafolder 'data.mat']);
handles.datafolder=datafolder;
p=p_original;
% load p parameters from data.mat in datafolder
set(handles.checkbox1,'Value',p.savetiff);
set(handles.checkbox2,'Value',p.savefinalresult);
set(handles.checkbox9,'Value',p.rawImage16);
set(handles.checkbox8,'Value',p.closeallfigure);
set(handles.popupmenu3,'Value',p.imageformat);
set(handles.edit6,'String',num2str(p.DAPI_chal));
set(handles.edit7,'String',num2str(p.FISH_chal));
if strcmp(p.gau_type,'gaussian')==1
    set(handles.popupmenu1,'Value',1);
else
    set(handles.popupmenu1,'Value',1);
end
set(handles.edit3,'String',num2str(p.gau_hsize));
set(handles.edit4,'String',num2str(p.gau_sigma));
set(handles.edit8,'String',num2str(p.init_stack));
set(handles.edit9,'String',num2str(p.final_stack));
set(handles.edit2,'String',num2str(p.back_th_num));
set(handles.edit5,'String',num2str(p.back_back_erode));
set(handles.edit10,'String',num2str(p.back_back_open));
set(handles.edit11,'String',num2str(p.back_nuc_open));
set(handles.edit12,'String',num2str(p.back_nuc_peakthresh));
set(handles.edit13,'String',num2str(p.candi_start));
set(handles.edit14,'String',num2str(p.candi_end));
set(handles.edit15,'String',num2str(p.nearby_radius_max));
set(handles.edit16,'String',num2str(p.remove_bou_back));
set(handles.edit17,'String',num2str(p.dist_radius_ratio));


save([datafolder 'data.mat'],'p','-append');
guidata(hObject, handles);


end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
load(['/functions/' 'root.mat']);
load([datafolder 'data.mat']);
handles.datafolder=datafolder;
p2=p;
p2.savetiff=get(handles.checkbox1,'Value');
 
p=p2;
save([datafolder 'data.mat'],'p','-append');
guidata(hObject, handles);
end

% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2
load(['/functions/' 'root.mat']);
load([datafolder 'data.mat']);
handles.datafolder=datafolder;
p2=p;
p2.savefinalresult=get(handles.checkbox2,'Value');
 
p=p2;
save([datafolder 'data.mat'],'p','-append');
guidata(hObject, handles);

end

function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
load(['/functions/' 'root.mat']);
load([datafolder 'data.mat']);
handles.datafolder=datafolder;
p2=p;
p2.gau_hsize=str2num(get(handles.edit3,'String'));
 
p=p2;
save([datafolder 'data.mat'],'p','-append');
guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double
load(['/functions/' 'root.mat']);
load([datafolder 'data.mat']);
handles.datafolder=datafolder;
p2=p;
p2.gau_sigma=str2num(get(handles.edit4,'String'));
 
p=p2;
save([datafolder 'data.mat'],'p','-append');
guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
load(['/functions/' 'root.mat']);
load([datafolder 'data.mat']);
handles.datafolder=datafolder;
p2=p;
p2.gau_type=get(handles.popupmenu1,'String');

p=p2;
save([datafolder 'data.mat'],'p','-append');
guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2
end

% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3
end

% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4
end


function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
end

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double
end

% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5
end

% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox6
end

% --- Executes on button press in checkbox7.
function checkbox7_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox7
end


% --- Executes on button press in checkbox8.
function checkbox8_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox8
load(['/functions/' 'root.mat']);
load([datafolder 'data.mat']);
handles.datafolder=datafolder;
p2=p;
p2.closeallfigure=get(handles.checkbox8,'Value');
 
p=p2;
save([datafolder 'data.mat'],'p','-append');
guidata(hObject, handles);
end

% --- Executes on button press in checkbox9.
function checkbox9_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox9
load(['/functions/' 'root.mat']);
load([datafolder 'data.mat']);
handles.datafolder=datafolder;
p2=p;
p2.rawImage16=get(handles.checkbox9,'Value');
 
p=p2;
save([datafolder 'data.mat'],'p','-append');
guidata(hObject, handles);
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double
load(['/functions/' 'root.mat']);
load([datafolder 'data.mat']);
handles.datafolder=datafolder;
p2=p;
p2.init_stack=max(min(str2num(get(handles.edit8,'String')),size(tiff_images,2)),1);
 
p=p2;
save([datafolder 'data.mat'],'p','-append');
guidata(hObject, handles);
end

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
end


function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double
load(['/functions/' 'root.mat']);
load([datafolder 'data.mat']);
handles.datafolder=datafolder;
p2=p;
p2.final_stack=max(min(str2num(get(handles.edit9,'String')),size(tiff_images,2)),1);
 
p=p2;
save([datafolder 'data.mat'],'p','-append');
guidata(hObject, handles);
end

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
end

% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3
load(['/functions/' 'root.mat']);
load([datafolder 'data.mat']);
handles.datafolder=datafolder;
p2=p;
p2.lsm_type=get(handles.popupmenu3,'String');
 
p=p2;
save([datafolder 'data.mat'],'p','-append');
guidata(hObject, handles);

end

% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double
load(['/functions/' 'root.mat']);
load([datafolder 'data.mat']);
handles.datafolder=datafolder;
p2=p;
p2.DAPI_chal=str2num(get(handles.edit6,'String'));
 
p=p2;
save([datafolder 'data.mat'],'p','-append');
guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double
load(['/functions/' 'root.mat']);
load([datafolder 'data.mat']);
handles.datafolder=datafolder;
p2=p;
p2.FISH_chal=str2num(get(handles.edit7,'String'));
 
p=p2;
save([datafolder 'data.mat'],'p','-append');
guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function edit17_Callback(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit17 as text
%        str2double(get(hObject,'String')) returns contents of edit17 as a double
end

% --- Executes during object creation, after setting all properties.
function edit17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double
end

% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double
end

% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double
end

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
end


function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit16 as text
%        str2double(get(hObject,'String')) returns contents of edit16 as a double
end

% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double
end

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
end


function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double
end

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
end


function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double
end

% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

