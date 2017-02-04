function varargout = movietest(varargin)
% MOVIETEST MATLAB code for movietest.fig
%      MOVIETEST, by itself, creates a new MOVIETEST or raises the existing
%      singleton*.
%
%      H = MOVIETEST returns the handle to a new MOVIETEST or the handle to
%      the existing singleton*.
%
%      MOVIETEST('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MOVIETEST.M with the given input arguments.
%
%      MOVIETEST('Property','Value',...) creates a new MOVIETEST or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before movietest_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to movietest_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help movietest

% Last Modified by GUIDE v2.5 12-Jul-2016 12:57:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @movietest_OpeningFcn, ...
    'gui_OutputFcn',  @movietest_OutputFcn, ...
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


% --- Executes just before movietest is made visible.
function movietest_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to movietest (see VARARGIN)

% Choose default command line output for movietest
handles.output = hObject;

% register mouse click event
set(gcf, 'WindowButtonDownFcn', @getMousePositionOnImage);
set(gcf, 'WindowButtonMotionFcn', @displayMousePosition);
pan off;

% Update handles structure
guidata(hObject, handles);

% Clear screen and variables
clear all;
clc;

function displayMousePosition(src, event)
handles = guidata(src);
    cursorPoint = get(handles.fig_movieDisp, 'CurrentPoint');
    X = cursorPoint(1,1);
    Y = cursorPoint(1,2);
    xLimits = get(handles.fig_movieDisp, 'xlim');
    yLimits = get(handles.fig_movieDisp, 'ylim');
    
    if (X> min(xLimits) && X < max(xLimits) && Y > min(yLimits) && Y < max(yLimits))
    
    set(handles.txt_mouseX,'String',['X: ' num2str(X)]);
    set(handles.txt_mouseY,'String',['Y: ' num2str(Y)]);
    end


function getMousePositionOnImage(src, event)
global mouseData numOfCurrentFrame algorithmFlag numOfTotalAnts
handles = guidata(src);

if algorithmFlag == 1
else
    
    cursorPoint = get(handles.fig_movieDisp, 'CurrentPoint');
    X = cursorPoint(1,1);
    Y = cursorPoint(1,2);
    
    set(handles.txt_mouseX,'String',['X: ' num2str(X)]);
    set(handles.txt_mouseY,'String',['Y: ' num2str(Y)]);
    
    xLimits = get(handles.fig_movieDisp, 'xlim');
    yLimits = get(handles.fig_movieDisp, 'ylim');
    
    clickType = get(src,'SelectionType');
    
    if (X> min(xLimits) && X < max(xLimits) && Y > min(yLimits) && Y < max(yLimits))
        antNum = str2double(get(handles.txt_antNum,'String'));
        if antNum <= numOfTotalAnts
            switch clickType
                case 'normal'
                    mouseData.coord(numOfCurrentFrame, antNum,:) = [X, Y, 1];
                    prev_list = get(handles.list_antNum,'String');
                    if ~ismember(num2str(antNum),prev_list)
                        new_list = sort([prev_list; num2str(antNum)]);
                        set(handles.list_antNum,'string',new_list);
                    end
                case 'alt'
                    interaction = get(handles.txt_interaction,'String');
                    [strAnts, strType] = strsplit(interaction,{'<','>','-'});
                    firstAnt = str2double(strAnts{1});
                    secondAnt = str2double(strAnts{2});
                    
                    if ~isnan(firstAnt) || ~isnan(secondAnt)
                        switch strType{1}
                            case '>'
                                interactionType = 1;
                            case '<'
                                interactionType = 2;
                            case '-'
                                interactionType = 3;
                            otherwise
                                interactionType = 0;
                        end
                        if firstAnt < 1 || firstAnt > numOfTotalAnts || ...
                                secondAnt < 1 || secondAnt > numOfTotalAnts || ...
                                firstAnt == secondAnt || interactionType == 0
                        else
                            dimension = size(mouseData.interaction);
                            duplicate = 0;
                            if numOfCurrentFrame > dimension(1)
                                mouseData.interaction(numOfCurrentFrame,end+1,:) = ...
                                    [X,Y,firstAnt,secondAnt,interactionType];
                            else
                                for i = 1:dimension(2)
                                    if mouseData.interaction(numOfCurrentFrame,i,3) == firstAnt && ...
                                            mouseData.interaction(numOfCurrentFrame,i,4) == secondAnt
                                        mouseData.interaction(numOfCurrentFrame,i,:) = ...
                                            [X,Y,firstAnt,secondAnt,interactionType];
                                        duplicate = 1;
                                        break;
                                    elseif mouseData.interaction(numOfCurrentFrame,i,4) == firstAnt && ...
                                            mouseData.interaction(numOfCurrentFrame,i,3) == secondAnt
                                        mouseData.interaction(numOfCurrentFrame,i,:) = ...
                                            [X,Y,firstAnt,secondAnt,interactionType];
                                        duplicate = 1;
                                        break;
                                    else
                                    end
                                end
                                if ~duplicate
                                    mouseData.interaction(numOfCurrentFrame,end+1,:) = ...
                                        [X,Y,firstAnt,secondAnt,interactionType];
                                end
                            end
                        end
                    end
                    
                otherwise
            end
            
            updatePlot(src,event);
        end
    else
    end
end

function updatePlot(src,event)
global CMap mouseData numOfCurrentFrame numOfTotalAnts
handles = guidata(src);
dimension = size(mouseData.coord);
dimension2 = size(mouseData.interaction);
interCount = str2double(get(handles.txt_interpolation,'String'));

if numOfCurrentFrame <= dimension(1)
    delete(findobj(handles.fig_movieDisp, 'type', 'line'));
    delete(findobj(handles.fig_movieDisp, 'type', 'text'));
    for i = 1:numOfTotalAnts
        if mouseData.coord(numOfCurrentFrame,i,3) == 1 || ...
              mouseData.coord(numOfCurrentFrame,i,3) == 0  
            if numOfCurrentFrame > 2
                for j = numOfCurrentFrame-1:-1:1
                    if mouseData.coord(j,i,3) == 1
                        if numOfCurrentFrame - j > 1
                            xEstimate = linspace(mouseData.coord(j,i,1),mouseData.coord(numOfCurrentFrame,i,1),numOfCurrentFrame-j+1);
                            xEstimate = xEstimate(2:end-1)';
                            if mouseData.coord(j,i,1) == mouseData.coord(numOfCurrentFrame,i,1)
                                x = [mouseData.coord(j,i,1) mouseData.coord(numOfCurrentFrame,i,1)+0.01];
                            else
                                x = [mouseData.coord(j,i,1) mouseData.coord(numOfCurrentFrame,i,1)];
                            end
                            if mouseData.coord(j,i,2) == mouseData.coord(numOfCurrentFrame,i,2)
                                y = [mouseData.coord(j,i,2) mouseData.coord(numOfCurrentFrame,i,2)+0.01];
                            else
                                y = [mouseData.coord(j,i,2) mouseData.coord(numOfCurrentFrame,i,2)];
                            end
                            yEstimate = interp1(x,y,xEstimate);
                            mouseData.coord(j+1:numOfCurrentFrame-1,i,:) = [xEstimate, yEstimate,ones(numOfCurrentFrame - j - 1,1)*2];
                        end
                        break;
                    end
                end
            end
            if numOfCurrentFrame < dimension(1) - 1
                for j = numOfCurrentFrame+1:dimension(1)
                    if mouseData.coord(j,i,3) == 1
                        if j - numOfCurrentFrame > 1
                            xEstimate = linspace(mouseData.coord(numOfCurrentFrame,i,1),mouseData.coord(j,i,1),j - numOfCurrentFrame+1);
                            xEstimate = xEstimate(2:end-1)';
                            if mouseData.coord(numOfCurrentFrame,i,1) == mouseData.coord(j,i,1)
                                x = [mouseData.coord(numOfCurrentFrame,i,1)+0.01 mouseData.coord(j,i,1)];
                            else
                                x = [mouseData.coord(numOfCurrentFrame,i,1) mouseData.coord(j,i,1)];
                            end
                            if mouseData.coord(numOfCurrentFrame,i,2) == mouseData.coord(j,i,2)
                                y = [mouseData.coord(numOfCurrentFrame,i,2)+0.01 mouseData.coord(j,i,2)];
                            else
                                y = [mouseData.coord(numOfCurrentFrame,i,2) mouseData.coord(j,i,2)];
                            end
                            yEstimate = interp1(x,y,xEstimate);
                            mouseData.coord(numOfCurrentFrame+1:j-1,i,:) = [xEstimate, yEstimate,ones(j - numOfCurrentFrame - 1,1)*2];
                        end
                        break;
                    end
                end
            end
        end
        
        if mouseData.coord(numOfCurrentFrame,i,3) == 0
            hold(handles.fig_movieDisp,'on');
            plot(handles.fig_movieDisp, mouseData.coord(numOfCurrentFrame,i,1),...
                mouseData.coord(numOfCurrentFrame,i,2),'^','Color',CMap(i,:),...
                'MarkerSize', 4);
            text(mouseData.coord(numOfCurrentFrame,i,1),...
                mouseData.coord(numOfCurrentFrame,i,2)+3,num2str(i),...
                'HorizontalAlignment','center','VerticalAlignment',...
                'top','Color',CMap(i,:));
            hold(handles.fig_movieDisp,'off');
        end
        
        if mouseData.coord(numOfCurrentFrame,i,3) == 1
            hold(handles.fig_movieDisp,'on');
            plot(handles.fig_movieDisp, mouseData.coord(numOfCurrentFrame,i,1),...
                mouseData.coord(numOfCurrentFrame,i,2),'o','Color',CMap(i,:),...
                'MarkerSize', 4);
            text(mouseData.coord(numOfCurrentFrame,i,1),...
                mouseData.coord(numOfCurrentFrame,i,2)+3,num2str(i),...
                'HorizontalAlignment','center','VerticalAlignment',...
                'top','Color',CMap(i,:));
            hold(handles.fig_movieDisp,'off');
        end
        if mouseData.coord(numOfCurrentFrame,i,3) == 2
            hold(handles.fig_movieDisp,'on');
            plot(handles.fig_movieDisp, mouseData.coord(numOfCurrentFrame,i,1),...
                mouseData.coord(numOfCurrentFrame,i,2),'*','Color',CMap(i,:),...
                'MarkerSize',4);
            text(mouseData.coord(numOfCurrentFrame,i,1),...
                mouseData.coord(numOfCurrentFrame,i,2)+3,num2str(i),...
                'HorizontalAlignment','center','VerticalAlignment',...
                'top','Color',CMap(i,:));
            hold(handles.fig_movieDisp,'off');
        end
%         if interCount > numOfCurrentFrame
%             index = find(mouseData.coord(1:numOfCurrentFrame,i,3) == 1|...
%                 mouseData.coord(1:numOfCurrentFrame,i,3) == 2,interCount,'Last');
%             hold(handles.fig_movieDisp,'on');
%             plot(handles.fig_movieDisp, mouseData.coord(index,i,1),...
%                 mouseData.coord(index,i,2),'-','Color',CMap(i,:));
%             hold(handles.fig_movieDisp,'off');
%         else
            index = find(mouseData.coord(1:numOfCurrentFrame,i,3) == 1|...
                mouseData.coord(1:numOfCurrentFrame,i,3) == 2,interCount,'Last');
            hold(handles.fig_movieDisp,'on');
            plot(handles.fig_movieDisp, mouseData.coord(index,i,1),...
                mouseData.coord(index,i,2),'-','Color',CMap(i,:));
            hold(handles.fig_movieDisp,'off');
%         end
    end
    
    if numOfCurrentFrame <= dimension2(1)
        for i = 2:dimension2(2)
            switch mouseData.interaction(numOfCurrentFrame,i,5)
                case 1
                    strType = '>>';
                case 2
                    strType = '<<';
                case 3
                    strType = '<->';
                otherwise
                    strType = '##';
            end
            hold(handles.fig_movieDisp,'on');
            if ~strcmp(strType, '##')
                plot(handles.fig_movieDisp, mouseData.interaction(numOfCurrentFrame,i,1),...
                    mouseData.interaction(numOfCurrentFrame,i,2),'x','Color','black',...
                    'MarkerSize', 4);
                text(mouseData.interaction(numOfCurrentFrame,i,1),...
                    mouseData.interaction(numOfCurrentFrame,i,2)-3,...
                    [num2str(mouseData.interaction(numOfCurrentFrame,i,3)), ...
                    strType, num2str(mouseData.interaction(numOfCurrentFrame,i,4))],...
                    'HorizontalAlignment','center',...
                    'VerticalAlignment','bottom',...
                    'Color','black');
            end
            hold(handles.fig_movieDisp,'off');
        end
    end
end


% --- Outputs from this function are returned to the command line.
function varargout = movietest_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in btn_selFile.
function btn_selFile_Callback(hObject, eventdata, handles)
% hObject    handle to btn_selFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global numOfTotalFrames numOfCurrentFrame fileList PathName ...
    mouseData CMap algorithmFlag selectedList backgroundFlag ...
    regionSelectFlag brightD BW numOfTotalAnts area antSize
[FileName,PathName,FilterIndex] = uigetfile({'*.jpg';'*.tif';'*.png'},'Choose a movie file');
if isequal(FileName,0)
    return
end
% get all files with the same extension
[prefix,junk,ext]=fileparts(FileName);
fileList=dir([PathName ['*' ext]]);

% get the total number of available files
numOfTotalFrames = numel(fileList);
set(handles.txt_totalFrameNum, 'String', num2str(numOfTotalFrames));

% get number of total ants in the frame
numOfTotalAnts = round(str2double(get(handles.txt_maxAntNum,'String')));

% set color map depending on the number of ants
CMap=colormap(hsv(numOfTotalAnts));

% set slider value to 1
set(handles.sli_navMovie, ...
    'value',1, ...
    'max',numOfTotalFrames, ...
    'min',1, ...
    'sliderstep',[1, 10]/numOfTotalFrames);

% set current
numOfCurrentFrame = 1;
antSize = str2double(get(handles.txt_antSize,'String'));

% update current frame number
set(handles.txt_currentFrameNum, 'String', num2str(numOfCurrentFrame));

% show the first frame image
axes(handles.fig_movieDisp);
imshow([PathName fileList(numOfCurrentFrame).name]);

% get list of ants
selectedList = get(handles.list_antNum,{'string','value'});

% initialize mouseData structure
mouseData.coord(1,1,:) = [0, 0, 0];
mouseData.interaction(1,1,:) = [0, 0, 0, 0, 0];

% initialize algorithm flag
algorithmFlag = 0;
backgroundFlag = 0;
regionSelectFlag = 0;

if exist([PathName 'savedFiles/bright_image.jpg']) == 2
    brightD = imread([PathName 'savedFiles/bright_image.jpg']);
    backgroundFlag = 1;
end
if exist([PathName 'savedFiles/BW.jpg']) == 2
    BW = imread([PathName 'savedFiles/BW.jpg']);
    area = bwarea(BW);
    regionSelectFlag = 1;
end
btn_trackingAlgorithm_Callback(hObject,eventdata,handles);

% --- Executes on slider movement.
function sli_navMovie_Callback(hObject, eventdata, handles)
% hObject    handle to sli_navMovie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global numOfCurrentFrame fileList PathName
numOfCurrentFrame = round(get(hObject,'value'));
axes(handles.fig_movieDisp);
imshow([PathName fileList(numOfCurrentFrame).name]);
set(handles.txt_currentFrameNum,'string', num2str(numOfCurrentFrame));
btn_trackingAlgorithm_Callback(hObject,eventdata,handles);


% --- Executes during object creation, after setting all properties.
function sli_navMovie_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sli_navMovie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function txt_currentFrameNum_Callback(hObject, eventdata, handles)
% hObject    handle to txt_currentFrameNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_currentFrameNum as text
%        str2double(get(hObject,'String')) returns contents of txt_currentFrameNum as a double
global mov1 numOfTotalFrames numOfCurrentFrame fileList PathName
numOfCurrentFrame = str2double(get(hObject,'String'));
if numOfCurrentFrame < 1
    numOfCurrentFrame = 1;
    set(hObject,'String',num2str(numOfCurrentFrame));
elseif numOfCurrentFrame > numOfTotalFrames
    numOfCurrentFrame = numOfTotalFrames;
    set(hObject,'String',num2str(numOfCurrentFrame));
end
set(handles.sli_navMovie,'value',numOfCurrentFrame);
axes(handles.fig_movieDisp);
imshow([PathName fileList(numOfCurrentFrame).name]);
btn_trackingAlgorithm_Callback(hObject,eventdata,handles);

% --- Executes during object creation, after setting all properties.
function txt_currentFrameNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_currentFrameNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on key press with focus on mainWindow and none of its controls.
function mainWindow_KeyPressFcn(hObject, eventdata, handles)
global numOfTotalFrames numOfCurrentFrame fileList PathName

timeStep = str2double(get(handles.txt_timeStep,'String'));
switch eventdata.Key
    case 'd'
        newTime = numOfCurrentFrame + timeStep;
    case 'a'
        newTime = numOfCurrentFrame - timeStep;
    case 's'
        antNum = str2double(get(handles.txt_antNum,'String'));
        if antNum > 1
            antNum = antNum - 1;
        end
        set(handles.txt_antNum, 'String', num2str(antNum));
        return
    case 'w'
        antNum = str2double(get(handles.txt_antNum,'String'));
        maxNum = str2double(get(handles.txt_maxAntNum,'String'));
        if antNum < maxNum
            antNum = antNum + 1;
        end
        set(handles.txt_antNum, 'String', num2str(antNum));
        return
    otherwise
        return
end
if newTime < 1
    newTime = 1;
elseif newTime > numOfTotalFrames
    newTime = numOfTotalFrames;
end
numOfCurrentFrame = newTime;
set(handles.txt_currentFrameNum, 'string', num2str(numOfCurrentFrame));
set(handles.sli_navMovie,'value',numOfCurrentFrame);
axes(handles.fig_movieDisp);
imshow([PathName fileList(numOfCurrentFrame).name]);
btn_trackingAlgorithm_Callback(hObject,eventdata,handles);

function calculateTrajectory(hObject, eventdata)
global mouseData numOfCurrentFrame BW brightD fileList PathName...
    selectedList numOfTotalAnts antSize
handles = guidata(hObject);
threshold = str2double(get(handles.txt_threshold,'String'));
dimension = size(mouseData.coord);

% Calculate K Means

dOriginal = imread([PathName fileList(numOfCurrentFrame).name]);
if size(dOriginal,3)==3
    d = rgb2gray(dOriginal);
else
    d = dOriginal;
end
diffD = (brightD - d) .* BW;
diffD = diffD * 2;
diffD = medfilt2(diffD);
[ki,kj,~] = ind2sub(size(diffD), find(diffD > threshold));

% for calculating area occupied by ants and estimate the number of ants
diffDbw = imextendedmax(diffD,threshold);
detectedArea = bwarea(diffDbw);

%[ki,kj] = find(diffDbw);
%disp(detectedArea);

assignin('base','diffD',diffD);
assignin('base','d',d);
%assignin('base','ki',ki);
%assignin('base','kj',kj);
%assignin('base','antSize',antSize);
%numOfTotalAnts = round(length(ki)/antSize);
%set(handles.txt_maxAntNum,'String',num2str(numOfTotalAnts));

% First frame should not be evaluated
if numOfCurrentFrame == 1
    %if(length(ki) > numOfTotalAnts*2)
        [~, C] = kmeans([kj ki], numOfTotalAnts);
    %end
    mouseData.coord(numOfCurrentFrame, 1:numOfTotalAnts,:) = [C(:,1), C(:,2), ones(numOfTotalAnts,1)];
else
    

    if numOfCurrentFrame <= dimension(1) && ...
            sum(find(mouseData.coord(numOfCurrentFrame,:,3) == 1)) &&...
            dimension(2) >= numOfTotalAnts        
    else
        if dimension(2) < numOfTotalAnts
            [~, C] = kmeans([kj ki], numOfTotalAnts);
            mouseData.coord(numOfCurrentFrame, 1:numOfTotalAnts,:) = [C(:,1), C(:,2), ones(numOfTotalAnts,1)];
        else
            
            preDist = [];
            % search backward on mouseData.coord for user input coordinates
              searchIndex = min(dimension(1), numOfCurrentFrame);
            for i = 1:numOfTotalAnts
                indX = find(mouseData.coord(1:searchIndex,i,3) == 1, 1, 'Last');
                if ~isempty(indX)
                    indY = i;
                    preDist(i,1) = mouseData.coord(indX,indY,1);
                    preDist(i,2) = mouseData.coord(indX,indY,2);
                end
            end
            
            if(length(ki) > numOfTotalAnts*2)
                [~, C] = kmeans([kj ki], numOfTotalAnts,'start', preDist);
            end
            
            if ~exist('C','var')
                C = preDist;
            end
            
            [ind, distance] = knnsearch(C, preDist);
            u = unique(ind);
            [n, b] = histc(ind, u);
            singleInd = find(ismember(b, find(n == 1)));
            ind1 = find(n>1);
            
            if max(n) > 1
                for i = 1:numel(ind1)
                    ind2 = find(ind == u(ind1(i)));
                    [t, ind3] = min(distance(ind2));
                    singleInd = [singleInd; ind2(ind3)];
                end
            end
            
            missingInd = 1:numOfTotalAnts;
            missingInd(singleInd) = [];
            missingVal = 1:numOfTotalAnts;
            missingVal(u) = [];
            ind(missingInd) = missingVal;
            
            C = C(ind,:);
            oldInd = [];
            for i = 1:length(selectedList{2})
                if str2num(selectedList{1}{selectedList{2}(i)}) > 0
                    oldInd = [oldInd,str2num(selectedList{1}{selectedList{2}(i)})];
                end
            end
            C(oldInd,:) = preDist(oldInd,:);
            mouseData.coord(numOfCurrentFrame, oldInd,:) =  [C(oldInd,1), C(oldInd,2), zeros(length(oldInd),1)];
            newInd = setdiff(1:numOfTotalAnts,oldInd);
            for i = 1:length(newInd)
                if sqrt(sum((C(newInd(i),:) - preDist(newInd(i),:)).^2)) > 70
                    C(newInd(i),:) = preDist(newInd(i),:);
                end
            end
            mouseData.coord(numOfCurrentFrame, newInd,:) = [C(newInd,1), C(newInd,2), ones(length(newInd),1)];
            
        end
    end
end

function updateAntList(hObject, eventdata)
global numOfTotalAnts selectedList
handles = guidata(hObject);
new_list = {'0'};
for i = 1:numOfTotalAnts
    new_list = sort([new_list, num2str(i)]);
    %prev_list = get(handles.list_antNum,'String');
    %if ~ismember(num2str(i),prev_list)
    %new_list = sort([prev_list; num2str(i)]);
    set(handles.list_antNum,'String',new_list);
    %end
end
selectedList = get(handles.list_antNum,{'string','value'});

% --- Executes during object creation, after setting all properties.
function mainWindow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mainWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function txt_timeStep_Callback(hObject, eventdata, handles)
% hObject    handle to txt_timeStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_timeStep as text
%        str2double(get(hObject,'String')) returns contents of txt_timeStep as a double


% --- Executes during object creation, after setting all properties.
function txt_timeStep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_timeStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function txt_totalFrameNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_totalFrameNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function txt_antNum_Callback(hObject, eventdata, handles)
% hObject    handle to txt_antNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_antNum as text
%        str2double(get(hObject,'String')) returns contents of txt_antNum as a double


% --- Executes during object creation, after setting all properties.
function txt_antNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_antNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_maxAntNum_Callback(hObject, eventdata, handles)
% hObject    handle to txt_maxAntNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CMap numOfTotalAnts
% Hints: get(hObject,'String') returns contents of txt_maxAntNum as text
%        str2double(get(hObject,'String')) returns contents of txt_maxAntNum as a double

numOfTotalAnts = round(str2double(get(hObject,'String')));
CMap=colormap(hsv(numOfTotalAnts));
btn_trackingAlgorithm_Callback(hObject,eventdata,handles);

% --- Executes during object creation, after setting all properties.
function txt_maxAntNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_maxAntNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in list_antNum.
function list_antNum_Callback(hObject, eventdata, handles)
% hObject    handle to list_antNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns list_antNum contents as cell array
%        contents{get(hObject,'Value')} returns selected item from list_antNum
global selectedList
selectedList = get(hObject,{'string','value'});
% disp(selectedList{2});
% disp(selectedList{1});
% for i = 1:length(selectedList{2})
%     disp(selectedList{1}(selectedList{2}(i)))
% end

% --- Executes during object creation, after setting all properties.
function list_antNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to list_antNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function txt_interpolation_Callback(hObject, eventdata, handles)
% hObject    handle to txt_interpolation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_interpolation as text
%        str2double(get(hObject,'String')) returns contents of txt_interpolation as a double
updatePlot(hObject,eventdata);

% --- Executes during object creation, after setting all properties.
function txt_interpolation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_interpolation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_interaction_Callback(hObject, eventdata, handles)
% hObject    handle to txt_interaction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_interaction as text
%        str2double(get(hObject,'String')) returns contents of txt_interaction as a double


% --- Executes during object creation, after setting all properties.
function txt_interaction_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_interaction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_saveData.
function btn_saveData_Callback(hObject, eventdata, handles)
% hObject    handle to btn_saveData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global PathName mouseData
mouseData.coord(:,:,3) = 1;
save([PathName 'savedFiles/' 'antData.mat'], 'mouseData');

% --- Executes during object creation, after setting all properties.
function txt_progress_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_progress (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

function txt_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to txt_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_threshold as text
%        str2double(get(hObject,'String')) returns contents of txt_threshold as a double
btn_trackingAlgorithm_Callback(hObject,eventdata,handles);

% --- Executes during object creation, after setting all properties.
function txt_threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_bgCalculate.
function btn_bgCalculate_Callback(hObject, eventdata, handles)
% hObject    handle to btn_bgCalculate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global fileList PathName  brightD backgroundFlag area BW
if ~exist([PathName 'savedFiles'])
    mkdir(PathName,'savedFiles');
end
timelapseFlag = get(handles.radio_timeLapse,'Value');
proportionFlag = get(handles.radio_proportion,'Value');
if(timelapseFlag)
    darkV = VideoWriter([PathName 'savedFiles/timelapse.mp4'],'MPEG-4');
    darkV.FrameRate = 3;
    open(darkV);
end
h = waitbar(0,'Calulating Background'); 
for i = 1:numel(fileList)
    dOriginal = imread([PathName fileList(i).name]);
    if size(dOriginal,3)==3
        d = rgb2gray(dOriginal);
    else
        d = dOriginal;
    end
    
    if i==1
        brightD = d;
        if(timelapseFlag)
            darkD = d;
        end
        if(proportionFlag)
            subD = brightD - darkD;
            subD = subD .* BW;
            occupied = numel(find(subD >= 50));
            proportion = occupied / area;
        end
        % initialize brightD
    else
        a = find(d>brightD);
        brightD(a) = d(a);
        
        if(timelapseFlag)
            b = find(d<darkD);
            darkD(b) = d(b);
        end
        if(proportionFlag)
            subD = brightD - darkD;
            subD = subD .* BW;
            occupied = numel(find(subD >= 50));
            proportion = [proportion occupied / area];
        end
    end
    if mod(i,10) == 0
        writeVideo(darkV,darkD);
    end
    waitbar(i/numel(fileList),h,'Calulating Background'); 
end
close(h);
imwrite(brightD,[PathName 'savedFiles/bright_image.jpg']);
close(darkV);

backgroundFlag = 1;

% --- Executes on button press in btn_trackMask.
function btn_trackMask_Callback(hObject, eventdata, handles)
% hObject    handle to btn_trackMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global BW PathName fileList regionSelectFlag area algorithmFlag backgroundFlag
if ~exist([PathName 'savedFiles'])
    mkdir(PathName,'savedFiles');
end
dOriginal = imread([PathName fileList(1).name]);
if size(dOriginal,3)==3
    d = rgb2gray(dOriginal);
else
    d = dOriginal;
end
algorithmFlag = 1;
[BW, xi1, yi1] = roipoly(dOriginal);   % main region of interest
area = polyarea(xi1, yi1);    % calculate area
BW = uint8(BW);
imwrite(BW,[PathName 'savedFiles/BW.jpg']);
algorithmFlag = 0;
regionSelectFlag = 1;

% --- Executes on button press in btn_nontrackMask.
function btn_nontrackMask_Callback(hObject, eventdata, handles)
% hObject    handle to btn_nontrackMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global BW PathName fileList regionSelectFlag area algorithmFlag backgroundFlag
if regionSelectFlag == 0
else
    BW = imread([PathName 'savedFiles/BW.jpg']);
    dOriginal = imread([PathName fileList(1).name]);
    if size(dOriginal,3)==3
        d = rgb2gray(dOriginal);
    else
        d = dOriginal;
    end
    algorithmFlag = 1;
    [BWNot, xi1, yi1] = roipoly(dOriginal);   % main region of interest
    area = area - polyarea(xi1, yi1);    % calculate area
    BWNot = uint8(BWNot);
    BW = BW&(~BWNot);
    BW = uint8(BW);
    
    imwrite(BW,[PathName 'savedFiles/BW.jpg']);
    algorithmFlag = 0;
end


% --- Executes on button press in btn_trackingAlgorithm.
function btn_trackingAlgorithm_Callback(hObject, eventdata, handles)
% hObject    handle to btn_trackingAlgorithm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global backgroundFlag regionSelectFlag
if backgroundFlag == 1 && regionSelectFlag == 1
    updateAntList(hObject,eventdata);
    calculateTrajectory(hObject,eventdata);
    updatePlot(hObject,eventdata);
end


% --- Executes on button press in btn_createMovie.
function btn_createMovie_Callback(hObject, eventdata, handles)
% hObject    handle to btn_createMovie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global PathName numOfTotalFrames mouseData numOfTotalAnts CMap fileList
analyzedMovie = VideoWriter([PathName 'savedFiles/trackedMovie.mp4'],'MPEG-4');
analyzedMovie.FrameRate = 10;
open(analyzedMovie);
axes(handles.fig_movieDisp);
%h = waitbar(0,'Creating Movie');
for i = 1:numOfTotalFrames
    %waitbar(i/numOfTotalFrames,h,'Creating Movie');
    d = imread([PathName fileList(i).name]);
    imshow(d);
    hold(handles.fig_movieDisp,'on');
    for j = 1:numOfTotalAnts
        plot(handles.fig_movieDisp, mouseData.coord(i,j,1),...
            mouseData.coord(i,j,2),'o','Color',CMap(j,:),...
            'MarkerSize', 4);
        text(mouseData.coord(i,j,1),...
            mouseData.coord(i,j,2)+3,num2str(j),...
            'HorizontalAlignment','center','VerticalAlignment',...
            'top','Color',CMap(j,:));
    end
    hold(handles.fig_movieDisp,'off');
    %drawnow;
    currframe = getframe;
    writeVideo(analyzedMovie,currframe); 
end
%close(h);

% --- Executes on button press in radio_timeLapse.
function radio_timeLapse_Callback(hObject, eventdata, handles)
% hObject    handle to radio_timeLapse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_timeLapse


% --- Executes on button press in radio_proportion.
function radio_proportion_Callback(hObject, eventdata, handles)
% hObject    handle to radio_proportion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_proportion


% --- Executes on button press in btn_loadData.
function btn_loadData_Callback(hObject, eventdata, handles)
% hObject    handle to btn_loadData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global numOfTotalAnts
[baseFileName, basePathName] = uigetfile('*.mat', 'Select a mat file');
if isequal(baseFileName,0)
    return
end
load([basePathName baseFileName]);
dimension = size(mouseData.coord);
numOfTotalAnts = dimension(2);
set(handles.txt_maxAntNum,'String',num2str(numOfTotalAnts));
updateAntList(hObject,eventdata);



function txt_antSize_Callback(hObject, eventdata, handles)
% hObject    handle to txt_antSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_antSize as text
%        str2double(get(hObject,'String')) returns contents of txt_antSize as a double
global antSize
antSize = str2num(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function txt_antSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_antSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function txt_mouseX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_mouseX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function txt_mouseY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_mouseY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
