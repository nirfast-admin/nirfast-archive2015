function varargout = gui_set_values(varargin)
% GUI_SET_VALUES M-file for gui_set_values.fig
%      GUI_SET_VALUES, by itself, creates a new GUI_SET_VALUES or raises the existing
%      singleton*.
%
%      H = GUI_SET_VALUES returns the handle to a new GUI_SET_VALUES or the handle to
%      the existing singleton*.
%
%      GUI_SET_VALUES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_SET_VALUES.M with the given input arguments.
%
%      GUI_SET_VALUES('Property','Value',...) creates a new GUI_SET_VALUES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_set_values_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_set_values_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_set_values

% Last Modified by GUIDE v2.5 10-Oct-2011 08:56:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_set_values_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_set_values_OutputFcn, ...
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


% --- Executes just before gui_set_values is made visible.
function gui_set_values_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_set_values (see VARARGIN)

% Choose default command line output for gui_set_values
handles.output = hObject;
set(hObject,'Name','Set Mesh Values');
if(nargin > 3)
    for index = 1:2:(nargin-3),
        if nargin-3==index, break, end
        switch lower(varargin{index})
         case 'type'
         handles.type = varargin{index+1};
        end
    end
end

% Modify GUI based on type
if ~strcmp(handles.type,'fluor') && ~strcmp(handles.type,'fluor_bem')
    set(handles.muax,'Enable','off');
    set(handles.musx,'Enable','off');
    set(handles.muam,'Enable','off');
    set(handles.musm,'Enable','off');
    set(handles.muaf,'Enable','off');
    set(handles.eta,'Enable','off');
    set(handles.tau,'Enable','off');
end
if ~strcmp(handles.type,'spec') && ~strcmp(handles.type,'spec_bem')
    set(handles.sa,'Enable','off');
    set(handles.sp,'Enable','off');
    set(handles.chromophores,'Enable','off');
else
    set(handles.ri,'Enable','off');
end
if strcmp(handles.type,'stnd') == 0 && strcmp(handles.type,'stnd_spn') == 0 ...
    && strcmp(handles.type,'stnd_bem') == 0
    set(handles.mua,'Enable','off');
    set(handles.mus,'Enable','off');
end
if strcmp(handles.type,'stnd_spn') == 0
    set(handles.g,'Enable','off');
end

% find workspace variables
vars = evalin('base','whos;');
varnames = {};
[nv,junk] = size(vars);
nflag = 1;
for i=1:1:nv
    flag = evalin('base',strcat('isfield(',vars(i).name,',''type'')'));
    if flag && strcmp(handles.type,evalin('base',strcat(vars(i).name,'.type')))
        varnames{nflag} = vars(i).name;
        nflag = nflag + 1;
    end
end
if ~isempty(varnames)
    set(handles.variables_mesh,'String',varnames);
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_set_values wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_set_values_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function mesh_Callback(hObject, eventdata, handles)
% hObject    handle to mesh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mesh as text
%        str2double(get(hObject,'String')) returns contents of mesh as a double


% --- Executes during object creation, after setting all properties.
function mesh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mesh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse_mesh.
function browse_mesh_Callback(hObject, eventdata, handles)
% hObject    handle to browse_mesh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fn,pn] = myuigetfile('*.node','Input Mesh');
if fn == 0
    return;
end
temp = [pn fn];
set(handles.mesh,'String',temp(1:end-5));

guidata(hObject, handles);



function savemeshto_Callback(hObject, eventdata, handles)
% hObject    handle to savemeshto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of savemeshto as text
%        str2double(get(hObject,'String')) returns contents of savemeshto as a double


% --- Executes during object creation, after setting all properties.
function savemeshto_CreateFcn(hObject, eventdata, handles)
% hObject    handle to savemeshto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse_savemeshto.
function browse_savemeshto_Callback(hObject, eventdata, handles)
% hObject    handle to browse_savemeshto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fn,pn] = myuiputfile('','Save Mesh To');
if fn == 0
    return;
end
set(handles.savemeshto,'String',[pn fn]);

guidata(hObject, handles);



function region_Callback(hObject, eventdata, handles)
% hObject    handle to region (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of region as text
%        str2double(get(hObject,'String')) returns contents of region as a double


% --- Executes during object creation, after setting all properties.
function region_CreateFcn(hObject, eventdata, handles)
% hObject    handle to region (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mua_Callback(hObject, eventdata, handles)
% hObject    handle to mua (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mua as text
%        str2double(get(hObject,'String')) returns contents of mua as a double


% --- Executes during object creation, after setting all properties.
function mua_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mua (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mus_Callback(hObject, eventdata, handles)
% hObject    handle to mus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mus as text
%        str2double(get(hObject,'String')) returns contents of mus as a double


% --- Executes during object creation, after setting all properties.
function mus_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ri_Callback(hObject, eventdata, handles)
% hObject    handle to ri (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ri as text
%        str2double(get(hObject,'String')) returns contents of ri as a double


% --- Executes during object creation, after setting all properties.
function ri_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ri (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sa_Callback(hObject, eventdata, handles)
% hObject    handle to sa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sa as text
%        str2double(get(hObject,'String')) returns contents of sa as a double


% --- Executes during object creation, after setting all properties.
function sa_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sp_Callback(hObject, eventdata, handles)
% hObject    handle to sp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sp as text
%        str2double(get(hObject,'String')) returns contents of sp as a double


% --- Executes during object creation, after setting all properties.
function sp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function chromophores_Callback(hObject, eventdata, handles)
% hObject    handle to chromophores (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of chromophores as text
%        str2double(get(hObject,'String')) returns contents of chromophores as a double


% --- Executes during object creation, after setting all properties.
function chromophores_CreateFcn(hObject, eventdata, handles)
% hObject    handle to chromophores (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function muax_Callback(hObject, eventdata, handles)
% hObject    handle to muax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of muax as text
%        str2double(get(hObject,'String')) returns contents of muax as a double


% --- Executes during object creation, after setting all properties.
function muax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to muax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function musx_Callback(hObject, eventdata, handles)
% hObject    handle to musx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of musx as text
%        str2double(get(hObject,'String')) returns contents of musx as a double


% --- Executes during object creation, after setting all properties.
function musx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to musx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function muam_Callback(hObject, eventdata, handles)
% hObject    handle to muam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of muam as text
%        str2double(get(hObject,'String')) returns contents of muam as a double


% --- Executes during object creation, after setting all properties.
function muam_CreateFcn(hObject, eventdata, handles)
% hObject    handle to muam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function musm_Callback(hObject, eventdata, handles)
% hObject    handle to musm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of musm as text
%        str2double(get(hObject,'String')) returns contents of musm as a double


% --- Executes during object creation, after setting all properties.
function musm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to musm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function muaf_Callback(hObject, eventdata, handles)
% hObject    handle to muaf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of muaf as text
%        str2double(get(hObject,'String')) returns contents of muaf as a double


% --- Executes during object creation, after setting all properties.
function muaf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to muaf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function eta_Callback(hObject, eventdata, handles)
% hObject    handle to eta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eta as text
%        str2double(get(hObject,'String')) returns contents of eta as a double


% --- Executes during object creation, after setting all properties.
function eta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tau_Callback(hObject, eventdata, handles)
% hObject    handle to tau (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tau as text
%        str2double(get(hObject,'String')) returns contents of tau as a double


% --- Executes during object creation, after setting all properties.
function tau_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tau (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in done.
function done_Callback(hObject, eventdata, handles)
% hObject    handle to done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mainGUIhandle = nirfast;
mainGUIdata  = guidata(mainGUIhandle);
content = get(mainGUIdata.script,'String');
batch = get(mainGUIdata.batch_mode,'Value');

varname = mygenvarname(get(handles.mesh,'String'));

if get(handles.mua,'String')
    content{end+1} = strcat('val.mua=',get(handles.mua,'String'),';');
    if ~batch
        evalin('base',content{end});
    end
end
if get(handles.mus,'String')
    content{end+1} = strcat('val.mus=',get(handles.mus,'String'),';');
    if ~batch
        evalin('base',content{end});
    end
end
if get(handles.g,'String')
    content{end+1} = strcat('val.g=',get(handles.g,'String'),';');
    if ~batch
        evalin('base',content{end});
    end
end
if get(handles.ri,'String')
    content{end+1} = strcat('val.ri=',get(handles.ri,'String'),';');
    if ~batch
        evalin('base',content{end});
    end
end
if get(handles.muax,'String')
    content{end+1} = strcat('val.muax=',get(handles.muax,'String'),';');
    if ~batch
        evalin('base',content{end});
    end
end
if get(handles.musx,'String')
    content{end+1} = strcat('val.musx=',get(handles.musx,'String'),';');
    if ~batch
        evalin('base',content{end});
    end
end
if get(handles.muam,'String')
    content{end+1} = strcat('val.muam=',get(handles.muam,'String'),';');
    if ~batch
        evalin('base',content{end});
    end
end
if get(handles.musm,'String')
    content{end+1} = strcat('val.musm=',get(handles.musm,'String'),';');
    if ~batch
        evalin('base',content{end});
    end
end
if get(handles.muaf,'String')
    content{end+1} = strcat('val.muaf=',get(handles.muaf,'String'),';');
    if ~batch
        evalin('base',content{end});
    end
end
if get(handles.eta,'String')
    content{end+1} = strcat('val.eta=',get(handles.eta,'String'),';');
    if ~batch
        evalin('base',content{end});
    end
end
if get(handles.tau,'String')
    content{end+1} = strcat('val.tau=',get(handles.tau,'String'),';');
    if ~batch
        evalin('base',content{end});
    end
end
if get(handles.sa,'String')
    content{end+1} = strcat('val.sa=',get(handles.sa,'String'),';');
    if ~batch
        evalin('base',content{end});
    end
end
if get(handles.sp,'String')
    content{end+1} = strcat('val.sp=',get(handles.sp,'String'),';');
    if ~batch
        evalin('base',content{end});
    end
end
if get(handles.chromophores,'String')
    [token1,remain] = strtok(get(handles.chromophores,'String'));
    [token2,remain] = strtok(remain);
    while token2
        content{end+1} = strcat('val.',token1,'=',token2,';');
        if ~batch
            evalin('base',content{end});
        end
        [token1,remain] = strtok(remain);
        [token2,remain] = strtok(remain);
    end
end

% FOR ALL  TYPES

meshloc = get_pathloc(get(handles.mesh,'String'));

if  evalin('base',['ischar(' meshloc ')'])
    foo = evalin('base',['load_mesh(' meshloc ')']);
else
    foo = evalin('base',meshloc);
end
nregions = unique(foo.region(:));

foo = get(handles.region,'String');
if isempty(foo) || isnan(str2double(foo)) || ...
        int16(str2double(foo)) ~= str2double(foo) || ...
        int16(str2double(foo)) < min(nregions) || ...
        int16(str2double(foo)) > max(nregions)
    errordlg(' You need to specify a valid region number!',...
        'Invalid Region Number')
    error(' You need to specify a region number!')
else
    rid_ = num2str(int16(str2double(foo)));
end

content{end+1} = strcat(varname,' = set_mesh(',meshloc,',',...
    rid_,',val);');
if ~batch
    evalin('base',content{end});
end
if get(handles.savemeshto,'String')
    savemeshto = get(handles.savemeshto,'String');
    if ~canwrite(savemeshto)
        [junk fn ext1] = fileparts(savemeshto);
        savemeshto = [tempdir fn ext1];
        disp(['No write access, writing here instead: ' savemeshto]);
    end
    
    content{end+1} = strcat('save_mesh(',varname,',''',...
        savemeshto,''');');
    if ~batch
        evalin('base',content{end});
    end
end

set(mainGUIdata.script, 'String', content);
guidata(nirfast, mainGUIdata);
close(gui_set_values);


% --- Executes on selection change in variables_mesh.
function variables_mesh_Callback(hObject, eventdata, handles)
% hObject    handle to variables_mesh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns variables_mesh contents as cell array
%        contents{get(hObject,'Value')} returns selected item from variables_mesh
contents = get(hObject,'String');
set(handles.mesh,'String',contents{get(hObject,'Value')});

% display regions as point clouds of different colors
axes(handles.axes);
hold on
mesh = evalin('base',contents{get(hObject,'Value')});
regions = unique(mesh.region);
colors = {'r','g','b','c','m','y','k','w'};
for i=1:length(regions)
    ri = regions(i);
    if i > 8
        color = 'k';
    else
        color = colors{i};
    end
    if mesh.dimension == 2
        plot(mesh.nodes(mesh.region==ri,1),mesh.nodes(mesh.region==ri,2), ...
            [color 'o'],'MarkerSize',3);
    elseif mesh.dimension == 3
        plot3(mesh.nodes(mesh.region==ri,1),mesh.nodes(mesh.region==ri,2), ...
            mesh.nodes(mesh.region==ri,3),[color 'o'],'MarkerSize',3);
    end
end
legend(num2str(regions(1:length(regions))));
hold off

% prepopulate fields with averages
if isfield(mesh,'mua')
    set(handles.mua,'String',num2str(mean(mesh.mua)));
end
if isfield(mesh,'mus')
    set(handles.mus,'String',num2str(mean(mesh.mus)));
end
if isfield(mesh,'sa')
    set(handles.sa,'String',num2str(mean(mesh.sa)));
end
if isfield(mesh,'sp')
    set(handles.sp,'String',num2str(mean(mesh.sp)));
end
if isfield(mesh,'muax')
    set(handles.muax,'String',num2str(mean(mesh.muax)));
end
if isfield(mesh,'musx')
    set(handles.musx,'String',num2str(mean(mesh.musx)));
end
if isfield(mesh,'muam')
    set(handles.muam,'String',num2str(mean(mesh.muam)));
end
if isfield(mesh,'musm')
    set(handles.musm,'String',num2str(mean(mesh.musm)));
end
if isfield(mesh,'muaf')
    set(handles.muaf,'String',num2str(mean(mesh.muaf)));
end
if isfield(mesh,'eta')
    set(handles.eta,'String',num2str(mean(mesh.eta)));
end
if isfield(mesh,'tau')
    set(handles.tau,'String',num2str(mean(mesh.tau)));
end
if isfield(mesh,'chromscattlist') && isfield(mesh,'conc')
    str_tmp = '';
    all_sol = char(mesh.chromscattlist);
    for i=1:size(all_sol,1)-2
        str_tmp = [str_tmp deblank(all_sol(i,:)) ' ' num2str(mean(mesh.conc(:,i))) ' '];
    end
    set(handles.chromophores,'String',str_tmp);
end


% --- Executes during object creation, after setting all properties.
function variables_mesh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to variables_mesh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function g_Callback(hObject, eventdata, handles)
% hObject    handle to g (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of g as text
%        str2double(get(hObject,'String')) returns contents of g as a double


% --- Executes during object creation, after setting all properties.
function g_CreateFcn(hObject, eventdata, handles)
% hObject    handle to g (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
