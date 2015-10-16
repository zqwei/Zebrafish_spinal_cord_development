function varargout = Cluster_View_0(varargin)
% CLUSTER_VIEW_0 MATLAB code for Cluster_View_0.fig
%      CLUSTER_VIEW_0, by itself, creates a new CLUSTER_VIEW_0 or raises the existing
%      singleton*.
%
%      H = CLUSTER_VIEW_0 returns the handle to a new CLUSTER_VIEW_0 or the handle to
%      the existing singleton*.
%
%      CLUSTER_VIEW_0('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CLUSTER_VIEW_0.M with the given input arguments.
%
%      CLUSTER_VIEW_0('Property','Value',...) creates a new CLUSTER_VIEW_0 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Cluster_View_0_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Cluster_View_0_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Cluster_View_0

% Last Modified by GUIDE v2.5 16-Oct-2015 10:41:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Cluster_View_0_OpeningFcn, ...
                   'gui_OutputFcn',  @Cluster_View_0_OutputFcn, ...
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


% --- Executes just before Cluster_View_0 is made visible.
function Cluster_View_0_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Cluster_View_0 (see VARARGIN)

% Choose default command line output for Cluster_View_0
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
set(handles.axesLMatThres, 'ytick',[], 'xtick', []);
set(handles.axesTrace, 'ytick',[], 'xtick', []);
set(handles.axesLMat, 'ytick',[], 'xtick', []);
set(handles.axesCLocation, 'ytick',[], 'xtick', [], 'ztick', []);

% UIWAIT makes Cluster_View_0 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Cluster_View_0_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenuDataset.
function popupmenuDataset_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuDataset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuDataset contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuDataset
handles.pushbuttonDataset.Enable = 'on';


% --- Executes during object creation, after setting all properties.
function popupmenuDataset_CreateFcn(hObject, eventdata, handles) %#ok<DEFNU>
% hObject    handle to popupmenuDataset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

datasetsNames = {'Data_Dre_E1_HuCH2BGCaMP6f_0_20141006_041947_corrected';...
               'Data_Dre_E1_BTXinjHuCH2BGCaMP6f_TL_20140818_045650_corrected_signal';...
               'Data_Dre_E1_HuCGCaMP6f-mnx1TagRFP_0-1_20150410_032910.corrected'};
hObject.String = datasetsNames(1:2);
           


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over popupmenuDataset.
function popupmenuDataset_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuDataset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on button press in pushbuttonDataset.
function pushbuttonDataset_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonDataset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hObject.Enable = 'off';
fileName = handles.popupmenuDataset.String{handles.popupmenuDataset.Value};
load(['../../TempDat/', fileName, '.mat'], 'dff', 'tracks', 'timePoints');
load(['../../TempDat/', fileName, '_PSDPeakTime.mat'], 'PSDPeakTime');
load(['../../TempDat/', fileName, '_LMatTime.mat'], 'LMatTime');
data.dff            = dff;
data.timePoints     = timePoints;
data.numTimePoints  = length(timePoints(1:end-1));
data.tracks         = tracks;
data.PSDPeakTime    = PSDPeakTime;
data.LMatTime       = LMatTime;
data.currTime       = 1;
data.LMatThres      = handles.sliderTres.Value;
data.normalizedData = handles.checkboxNormalized.Value;

handles.output.UserData  = data;
handles.sliderTime.Value = 0;
handles.textTime.String  = ['Time: ' num2str(data.currTime)];
updateFactorList(handles.output.UserData, handles)
updatePlots(handles.output.UserData, handles.listboxCellType.Value, handles);


% --- Executes during object creation, after setting all properties.
function axesTrace_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axesTrace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axesTrace


% --- Executes on slider movement.
function sliderTime_Callback(hObject, eventdata, handles)
% hObject    handle to sliderTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

if isfield(handles, 'output')
    handles.output.UserData.currTime = ceil(hObject.Value*handles.output.UserData.numTimePoints);
    if handles.output.UserData.currTime<1; handles.output.UserData.currTime=1; end
    handles.textTime.String = ['Time: ' num2str(handles.output.UserData.currTime)];
    handles.listboxCellType.Value = 1;
    updateFactorList(handles.output.UserData, handles)
    updatePlots(handles.output.UserData, handles.listboxCellType.Value, handles);
end


% --- Executes during object creation, after setting all properties.
function sliderTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on selection change in listboxCellType.
function listboxCellType_Callback(hObject, eventdata, handles)
% hObject    handle to listboxCellType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listboxCellType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listboxCellType
updatePlots(handles.output.UserData, handles.listboxCellType.Value, handles);


% --- Executes during object creation, after setting all properties.
function listboxCellType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxCellType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

if isfield(handles, 'output')
    updateFactorList(handles.output.UserData, handles)
end


% --- Executes on slider movement.
function sliderTres_Callback(hObject, eventdata, handles)
% hObject    handle to sliderTres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.output.UserData.LMatThres = floor(hObject.Value*100)/100;
handles.textThres.String = ['Thres: ' num2str(handles.output.UserData.LMatThres)];
% disp(['Thres: ' num2str(handles.output.UserData.LMatThres)])
updatePlots(handles.output.UserData, handles.listboxCellType.Value, handles);


% --- Executes during object creation, after setting all properties.
function sliderTres_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderTres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
hObject.Value = 0.3;


% --- Executes on button press in checkboxNormalized.
function checkboxNormalized_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxNormalized (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxNormalized
handles.output.UserData.normalizedData = hObject.Value;
updatePlots(handles.output.UserData, handles.listboxCellType.Value, handles);


function updatePlots(UserData, cellType, handles)
dff                 = UserData.dff;
timePoints          = UserData.timePoints;
nTime               = UserData.currTime;
slicedDFF           = dff(:,timePoints(nTime)+1:timePoints(nTime+1)); 
numNeuron           = size(slicedDFF, 1);
tracks              = squeeze(mean(UserData.tracks(:, timePoints(nTime)+1:timePoints(nTime+1), :), 2));
PSDPeakTime         = UserData.PSDPeakTime(:, nTime);
LMat                = UserData.LMatTime{nTime};
Thres               = UserData.LMatThres;
normalizedData      = UserData.normalizedData;

cFactor             =  [0         0.4470    0.7410
                        0.8500    0.3250    0.0980
                        0.9290    0.6940    0.1250
                        0.4940    0.1840    0.5560
                        0.4660    0.6740    0.1880
                        0.3010    0.7450    0.9330
                        0.6350    0.0780    0.1840];

% Loading Matrix
axes(handles.axesLMat);
imagesc(LMat, 'Parent', handles.axesLMat);
colormap(handles.axesLMat, jet); 
caxis(handles.axesLMat, [-1 1])
axis xy
box off
set(handles.axesLMat, 'ytick', 1:numNeuron, 'xtick', []);

% Loading Matrix -- Threshold
axes(handles.axesLMatThres);
imagesc(LMat>Thres, 'Parent', handles.axesLMatThres); 
colormap(handles.axesLMatThres, gray(2)); 
caxis(handles.axesLMatThres, [0 1])
axis xy
set(handles.axesLMatThres, 'ytick',[], 'xtick', []);

% Neuron traces
cla(handles.axesTrace,'reset')
axes(handles.axesTrace);

numOrder           = 19;
lenWindow          = 511;


hold on;
for nNeuron         = 1:numNeuron
    nSlicedDFF      = slicedDFF(nNeuron, :);
    if ~normalizedData
        h = plot(nSlicedDFF + nNeuron, '-', 'linewid', 0.5);
        text(30, nNeuron+0.5, num2str(nNeuron));
    else
%         h = plot((nSlicedDFF - mean(nSlicedDFF))/std(nSlicedDFF) + nNeuron * 5, '-', 'linewid', 0.5);
        stdNSlicedDFF    = (nSlicedDFF - mean(nSlicedDFF))/std(nSlicedDFF);
        filtedNSlicedDFF = stdNSlicedDFF;
%         filtedNSlicedDFF = stdNSlicedDFF - sgolayfilt(stdNSlicedDFF, numOrder, lenWindow);
%         filtedNSlicedDFF = stdNSlicedDFF - medfilt1(stdNSlicedDFF);
        h = plot( filtedNSlicedDFF + nNeuron * 5, '-', 'linewid', 0.5);
%         plot( stdNSlicedDFF + nNeuron * 5, '--k', 'linewid', 0.5);
        text(30, nNeuron*5+2.5, num2str(nNeuron));
    end
    if isnan(PSDPeakTime(nNeuron))
        set(h, 'color', 'k')
    elseif PSDPeakTime(nNeuron)==0 %< 1/100 || PSDPeakTime(nNeuron) > 0.25 
        set(h, 'color', [0.5 0.5 0.5])
    else
        [~, nCellFactor] = max(LMat(nNeuron, :));
        set(h, 'color', cFactor(nCellFactor, :));
        if cellType > 1 
            if LMat(nNeuron, cellType-1)>Thres
                set(h, 'color', [cFactor(cellType-1, :) 1.0], 'linewid', 2.0);
            else
                set(h, 'color', [cFactor(nCellFactor, :) 0.5]);
            end
        end
    end
    
end
xlim([0 size(slicedDFF, 2)]);
if ~normalizedData
    ylim([-1 numNeuron+2])
else
    ylim([-1 numNeuron+1]*5)
end
grid on
% grid(handles.axesTrace, 'XGrid', 'on', 'YGrid', 'on')
set(handles.axesTrace, 'yticklabel',[], 'xticklabel', []);
set(handles.axesTrace, 'color','k');
hold off

% Neuron location

cla(handles.axesCLocation,'reset')
axes(handles.axesCLocation);
hold on;
for nNeuron         = 1:numNeuron
    if isnan(PSDPeakTime(nNeuron))
        scatter3(tracks(nNeuron,1), tracks(nNeuron,2), tracks(nNeuron,3), 8*10,...
            'o', 'linewid', 1.0, 'MarkerEdgeColor', 'k', 'MarkerFacecolor', 'w');
    elseif PSDPeakTime(nNeuron) == 0
        scatter3(tracks(nNeuron,1), tracks(nNeuron,2), tracks(nNeuron,3), 8*10, ...
            'o', 'linewid', 1.0, 'MarkerEdgeColor', [0.5 0.5 0.5], 'MarkerFacecolor', 'w');
    else
        [~, nCellFactor] = max(LMat(nNeuron, :));
        h = scatter3(tracks(nNeuron,1), tracks(nNeuron,2), tracks(nNeuron,3), ...
                (8 + 10*max(LMat(nNeuron, nCellFactor),0))*10,...
                'o', 'linewid', 1.0);
        set(h, 'MarkerFacecolor', cFactor(nCellFactor, :), 'MarkerEdgeColor', 'w');
        if cellType > 1 && LMat(nNeuron, cellType-1)>Thres
            set(h, 'MarkerEdgeColor', 'k');
        end        
    end
    text(tracks(nNeuron,1), tracks(nNeuron,2), tracks(nNeuron,3), num2str(nNeuron))
end
hold off
xlabel('x')
ylabel('y')
zlabel('z')
view(handles.axesCLocation, 45, 45)
    

function updateFactorList(UserData, handles)
currTime       = UserData.currTime;
% currPSDPeakTime= UserData.PSDPeakTime(:, currTime);
LMat           = UserData.LMatTime{currTime};
numList        = 1;

handles.listboxCellType.String = {'All neuron'};

% if sum(isnan(currPSDPeakTime))>0
%     numList    = numList + 1;
%     handles.listboxCellType.String{numList} = 'Non-active neuron';
% end
% 
% if sum(currPSDPeakTime == 0)>0
%     numList    = numList + 1;
%     handles.listboxCellType.String{numList} = 'Low frequency neuron';
% end

if size(LMat, 2) > 1
    for nFactor = 1:size(LMat, 2)
        numList    = numList + 1;
        handles.listboxCellType.String{numList} = ['Factor #' num2str(nFactor) ' neuron'];
    end
end
