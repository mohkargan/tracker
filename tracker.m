function varargout = tracker(varargin)
% TRACKER MATLAB code for tracker.fig
%      TRACKER, by itself, creates a new TRACKER or raises the existing
%      singleton*.
%
%      H = TRACKER returns the handle to a new TRACKER or the handle to
%      the existing singleton*.
%
%      TRACKER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRACKER.M with the given input arguments.
%
%      TRACKER('Property','Value',...) creates a new TRACKER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before tracker_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to tracker_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help tracker

% Last Modified by GUIDE v2.5 25-Dec-2017 15:42:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tracker_OpeningFcn, ...
                   'gui_OutputFcn',  @tracker_OutputFcn, ...
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

% --- Executes just before tracker is made visible.
function tracker_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tracker (see VARARGIN)

% Choose default command line output for tracker
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
handles.written_video = '';

set(handles.text_status,'String','Please set the settings!');

global video_finish;
video_finish = 0;
global replay_init;
replay_init = 0;
% Update handles structure
guidata(hObject,handles);

% UIWAIT makes tracker wait for user response (see UIRESUME)
% uiwait(handles.figure1);
end

% --- Outputs from this function are returned to the command line.
function varargout = tracker_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end

% --- Executes on button press in pushbutton_settings.
function pushbutton_settings_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_settings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[handles.padding, handles.lambda, handles.output_sigma_factor, handles.defaultROI, ...
handles.isIP, handles.opticalFlow, handles.kernel, handles.pathname, ...
handles.filename, handles.features, handles.cell_size] = settings;
if strcmp(handles.padding ,'') || strcmp(handles.lambda ,'') || strcmp(handles.output_sigma_factor ,'') || strcmp(handles.defaultROI ,'') ||...
 strcmp(handles.isIP ,'') || strcmp(handles.opticalFlow ,'') || strcmp(handles.kernel ,'') || strcmp(handles.pathname ,'') ||...
  strcmp(handles.filename ,'') || strcmp(handles.features ,'')  

    if ~strcmp(handles.kernel.type, 'linear') || strcmp(handles.cell_size ,'')       
    msgbox('Setting info is inaccurate, Please set again!','Info Retrieval Error!','error');
    clearAll(hObject,handles);
    return;
    end
end
set(handles.text_status,'String','Got the settings info!');
set(handles.edit_video_start, 'Enable','on');
set(handles.edit_video_end, 'Enable','on');
set(handles.checkbox_frames, 'Enable','on');
set(handles.checkbox_video, 'Enable','on');
set(handles.pushbutton_play, 'Enable','on');
fullfilename = fullfile(handles.pathname,handles.filename);
vidReader = VideoReader(fullfilename);
handles.vidReader = vidReader;
set(handles.text_duration,'String',num2str(vidReader.Duration));
set(handles.edit_video_end,'String',num2str(vidReader.Duration));
handles.video_duration = vidReader.Duration;
guidata(hObject,handles);
end

function edit_video_start_Callback(hObject, eventdata, handles)
% hObject    handle to edit_video_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_video_start as text
%        str2double(get(hObject,'String')) returns contents of edit_video_start as a double
end

% --- Executes during object creation, after setting all properties.
function edit_video_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_video_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit_video_end_Callback(hObject, eventdata, handles)
% hObject    handle to edit_video_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_video_end as text
%        str2double(get(hObject,'String')) returns contents of edit_video_end as a double
end

% --- Executes during object creation, after setting all properties.
function edit_video_end_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_video_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in checkbox_frames.
function checkbox_frames_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_frames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_frames
end

% --- Executes on button press in checkbox_video.
function checkbox_video_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_video (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_video
end

% --- Executes on button press in pushbutton_play.
function pushbutton_play_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_play (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

start_time = str2double(get(handles.edit_video_start,'String'));
video_end = str2double(get(handles.edit_video_end,'String'));
if isnan(start_time) || isnan(video_end)
    msgbox('Please enter values as numerics!', 'Error','error');
    return;
end
if (start_time < 0) || (video_end > handles.video_duration) || (start_time > video_end)
    if (video_end - handles.video_duration ) > 0.01
        msgbox('Please enter valid beginning and ending!', 'Error','error');
        return;
    end
end

set(handles.pushbutton_play,'Enable','off');
set(handles.pushbutton_finish,'Enable','on');
write_frames = get(handles.checkbox_frames,'Value');
write_video = get(handles.checkbox_video,'Value');
set(handles.checkbox_frames,'Enable','off');
set(handles.checkbox_video,'Enable','off');

global video_finish;
%% get info from handles
pathname = handles.pathname;
filename = handles.filename;
padding = handles.padding;
kernel = handles.kernel;
lambda = handles.lambda;
output_sigma_factor = handles.output_sigma_factor;
if ~strcmp(kernel.type,'linear')
    interp_factor = handles.kernel.interp_factor;
    cell_size = handles.cell_size;
else
    cell_size = 1;
    interp_factor = 0.02;
end
features = handles.features;
defaultroi = handles.defaultROI;
isIP = handles.isIP;
opticalFlow = handles.opticalFlow;
vidReader = handles.vidReader;

%% tracker in here, not function
    set(handles.text_status,'String','Processing the video!');
    %% find total frame count
    framecount = floor(vidReader.Duration * vidReader.FrameRate) + 1;
    
    %% set starting time
    if start_time > 0
        vidReader.CurrentTime = start_time;
    end

    
    %% set roi
    set(handles.axes1,'Visible','on');
    if defaultroi == 1
        roistruct = load('MVI_mask.mat','-mat');
        roi = roistruct.BW;
    else
        %kullanýcýnýn alaný sol ile seçip, birleþtirip en son alana sað
        %týklatýp "create mask" demesi gerekmektedir
        set(handles.text_status,'String','Processing the video: Draw ROI!');
        im = readFrame(vidReader);
        im = imresize(im,0.5);
        roi = roipoly(im);
    end
    set(handles.text_status,'String','Processing the video: ROI acquired!');
    %% set saving path according to folder of video, correlation, features and
    filename_noext = strtok(filename,'.'); %filename without extention
    if write_video
        video_path = strcat(filename_noext,'_',kernel.type,'_');
    end
    save_path = sprintf('%s%s%s',pathname,'frames_',filename_noext);
    
    save_path = strcat(save_path,'_',kernel.type,'_');
    
    if features.hog 
        save_path = strcat(save_path,'hog');
        if write_video
            video_path = strcat(video_path,'hog');
        end
    else
        save_path = strcat(save_path,'gray');
        if write_video
            video_path = strcat(video_path,'gray');
        end
    end
    
    if defaultroi == 1
        save_path = strcat(save_path,'_defroi'); %default roi
        if write_video
            video_path = strcat(video_path,'_defroi'); %default roi
        end
    else
        save_path = strcat(save_path,'_userroi'); % user selected roi
        if write_video
            video_path = strcat(video_path,'_userroi'); % user selected roi
        end
    end
    if start_time == 0
        if write_video
            video_path = strcat(video_path,'_beginto');
        end
        save_path = strcat(save_path,'_beginto');
    else
        if write_video
            video_path = sprintf('%s_%d%s',video_path,floor(start_time),'to');
        end
        save_path = sprintf('%s_%d%s',save_path,floor(start_time),'to');
    end
    if abs(video_end - handles.video_duration ) < 0.01
        if write_video
            video_path = strcat(video_path,'end.mp4');
        end
        save_path = strcat(save_path,'end\');
    else
        if write_video
            video_path = sprintf('%s%d%s',video_path,floor(video_end),'.mp4');
        end
        save_path = sprintf('%s%d%s',save_path,floor(video_end),'\');
    end
    if write_video
        video_path = fullfile(pathname,video_path);
    end
    mkdir(save_path);
    time = 0;  %to calculate FPS
    tic()
    %% optical flow ile otomatik seçme   
    if strcmp(opticalFlow.type,'Lucas-Kanade')
        opticFlow = opticalFlowLK('NoiseThreshold',opticalFlow.noiseThreshold);
    else %if strcmp(opticalFlow.type,'Horn-Schunck')
        opticFlow = opticalFlowHS('MaxIteration',opticalFlow.maxIteration,'Smoothness',opticalFlow.smoothness);
    end
    im1 = readFrame(vidReader);
    frameGray1 = rgb2gray(im1);
    im = readFrame(vidReader);
    
    if isIP == 1,
        whilecount = 0;
        while psnr(im,im1) > 40,
            im = readFrame(vidReader);
            whilecount = whilecount + 1;
        end
    end
    
    frameGray = rgb2gray(im);
    frameGray = imresize(frameGray, 0.5);
    frameGray1 = imresize(frameGray1, 0.5);

    [count,x,y,width,height] = ofmod(opticFlow,frameGray1, frameGray);
    %% bu kýsýmda baþlangýç için iç içe geçmeleri engelle ve roi uygula
    cx = floor(x+width/2);
    cy = floor(y+height/2);
    roicount = 1;
    inroi = int16.empty;
    for i = 1:count,
    	if (roi(cx(i),cy(i)))
        	inroi(roicount) = i;
            roicount = roicount + 1;
        end
    end
    count = size(inroi,2);
         
    if count ~= 0,
    	%% check bounding rectangles
        tempRects = zeros(count,4);
        for i = 1:count
            tempRects(i,:) = [x(inroi(i)); y(inroi(i)); width(inroi(i)); height(inroi(i))];
        end
        tempRects = checkRects(tempRects);
        
        new_x = tempRects(:,1);
        new_y = tempRects(:,2);
        new_w = tempRects(:,3);
        new_h = tempRects(:,4);
        count = size(new_x,1);
        
        target_sz= zeros(count,2);
        pos = zeros(count,2);
        for i=1:count
            target_sz(i,:) = [new_w(i), new_h(i)];
            pos(i,:) = [new_x(i), new_y(i)]+ floor(target_sz(i,:)/2);%merkez
        end
    else
        target_sz = 0;
        pos = 0;
    end
    %% without bounding rect control
%                    new_x = zeros(count,1); new_y = zeros(count,1); new_w = zeros(count,1); new_h = zeros(count,1);
%                     for i = 1:count
%                         new_x(i) = x(inroi(i));
%                         new_y(i) = y(inroi(i));
%                         new_w(i) = width(inroi(i));
%                         new_h(i) = height(inroi(i));
%                     end
        target_sz= zeros(count,2);
        pos = zeros(count,2);
        for i=1:count
            target_sz(i,:) = [new_w(i), new_h(i)];
            pos(i,:) = [new_x(i), new_y(i)]+ floor(target_sz(i,:)/2);%merkez
        end
    %% 
       
    if ~isempty(pos),
        window_sz = zeros(count,2);
        yf = cell(count,1);
        cos_window = cell(count,1);
        for i=1:count
            %window size, taking padding into account
            window_sz(i,:) = floor(target_sz(i) * (1 + padding));

            %create regression labels, gaussian shaped, with a bandwidth
            %proportional to target size
            output_sigma = sqrt(prod(target_sz(i,:))) * output_sigma_factor / cell_size;
            yf{i} = fft2(gaussian_shaped_labels(output_sigma, floor(window_sz(i,:) / cell_size)));

            %store pre-computed cosine window
            cos_window{i} = hann(size(yf{i},1)) * hann(size(yf{i},2))';		
        end
    end
	
    
    resize_image = 1;
    frame = 1;
    issame = 0;
    init = 1;
    of_interval = 15; %optical flow detection interval
    del_interval = 10; % interval to check if a object is in the same place for 

    if ~isempty(pos),
        positions = zeros(count,framecount, 2);  
        patch = cell(count,1);
        zf = cell(count,1);
        kzf = cell(count,1);
        model_xf = cell(count,1);
        kf = cell(count,1);
        xf = cell(count,1);
        alphaf = cell(count,1);
        model_alphaf = cell(count,1);
        poscount = zeros(count,1);
        peak_value = zeros(count,1);
    end
    fID = fopen(strcat(save_path,'fps.txt'),'w');
	while hasFrame(vidReader),

        if vidReader.CurrentTime > video_end || video_finish
            time = time + toc();
            break;
        end
        set(handles.text_status,'String',sprintf('Processing the video: %.1f seconds',vidReader.CurrentTime));
        %% 20 frame boyunca ayni pozisyondaysa trackerdan sil
        if count ~=0,
            tempcount = count;
            tempposcount = poscount;
            forcount = 1;
            temptarget_sz = target_sz;
            temppos = pos;
            tempwindow_sz = window_sz;
            for i = 1:count,
                if (poscount(i)>del_interval),
                    patch(i-(count-tempcount),:) = [];
                    zf(i-(count-tempcount),:) = [];
                    kzf(i-(count-tempcount),:) = [];
                    model_xf(i-(count-tempcount),:) = [];
                    kf(i-(count-tempcount),:) = [];
                    xf(i-(count-tempcount),:) = [];
                    alphaf(i-(count-tempcount),:) = [];
                    model_alphaf(i-(count-tempcount),:) = [];
                    yf(i-(count-tempcount),:) = [];
                    cos_window(i-(count-tempcount),:) = [];
                    tempcount = tempcount - 1;
                else
                    tempposcount(forcount) = poscount(i);
                    temptarget_sz(forcount,:) = target_sz(i,:);
                    temppos(forcount,:) = pos(i,:);
                    tempwindow_sz(forcount,:) = window_sz(i,:);
                    forcount = forcount + 1;
                end
            end
        
            poscount = tempposcount(1:tempcount);
            target_sz = temptarget_sz(1:tempcount,:);
            pos = temppos(1:tempcount,:);
            window_sz = tempwindow_sz(1:tempcount,:);
            count = tempcount;
        end
		%% load image
        if init  == 0,
            im = readFrame(vidReader);
        end
        %% optical flow buraya
        if mod(frame,of_interval) == 0,

            im = readFrame(vidReader);
            if size(im,3) > 1,
                im = rgb2gray(im);
            end
            if resize_image,
                im = imresize(im, 0.5);
            end
            if isIP == 1,
                while psnr(im,imtut) > 40;
                    im = readFrame(vidReader);
                    if size(im,3) > 1,
                        im = rgb2gray(im);
                    end
                    if resize_image,
                        im = imresize(im, 0.5);
                    end
                end
            end
            if strcmp(opticalFlow.type,'Lucas-Kanade')
                opticFlow = opticalFlowLK('NoiseThreshold',opticalFlow.noiseThreshold);
            elseif strcmp(opticalFlow.type,'Lucas-Kanade')
                opticFlow = opticalFlowHS('MaxIteration',opticalFlow.maxIteration,'Smoothness',opticalFlow.smoothness,0);
            end
            [count,x,y,width,height] = ofmod(opticFlow,imtut, im);

            %% check if is in roi ROI 
            cx = floor(x+width/2);
            cy = floor(y+height/2);
            roicount = 1;
            inroi = int16.empty;
            for i = 1:count,
                if (roi(cx(i),cy(i))),
                    inroi(roicount) = i;
                    roicount = roicount + 1;
                end
            end
            count = size(inroi,2);
           
            if count ~= 0,
                %% check bounding rectangles
                tempRects = zeros(count,4);
                for i = 1:count
                    tempRects(i,:) = [x(inroi(i)); y(inroi(i)); width(inroi(i)); height(inroi(i))];
                end
                tempRects = checkRects(tempRects);

                new_x = tempRects(:,1);
                new_y = tempRects(:,2);
                new_w = tempRects(:,3);
                new_h = tempRects(:,4);
                count = size(new_x,1);

                target_sz= zeros(count,2);
                pos = zeros(count,2);
                for i=1:count
                    target_sz(i,:) = [new_w(i), new_h(i)];
                    pos(i,:) = [new_x(i), new_y(i)]+ floor(target_sz(i,:)/2);%merkez
                end

%                    %% without bounding rect control
%                    new_x = zeros(count,1); new_y = zeros(count,1); new_w = zeros(count,1); new_h = zeros(count,1);
%                     for i = 1:count
%                         new_x(i) = x(inroi(i));
%                         new_y(i) = y(inroi(i));
%                         new_w(i) = width(inroi(i));
%                         new_h(i) = height(inroi(i));
%                     end
                    
                    target_sz= zeros(count,2);
                    pos = zeros(count,2);
                    for i=1:count
                        target_sz(i,:) = [new_w(i), new_h(i)];
                        pos(i,:) = [new_x(i), new_y(i)]+ floor(target_sz(i,:)/2);%merkez
                    end
                %%
                window_sz = zeros(count,2);
                yf = cell(count,1);
                cos_window = cell(count,1);
                for i=1:count
                    %window size, taking padding into account
                    window_sz(i,:) = floor(target_sz(i,:) * (1 + padding));

                    %create regression labels, gaussian shaped, with a bandwidth
                    %proportional to target size
                    output_sigma = sqrt(prod(target_sz(i,:))) * output_sigma_factor / cell_size;
                    yf{i} = fft2(gaussian_shaped_labels(output_sigma, floor(window_sz(i,:) / cell_size)));

                    %store pre-computed cosine window
                    cos_window{i} = hann(size(yf{i},1)) * hann(size(yf{i},2))';		
                end

                positions = zeros(count,framecount, 2);  

                patch = cell(count,1);
                zf = cell(count,1);
                kzf = cell(count,1);
                model_xf = cell(count,1);
                kf = cell(count,1);
                xf = cell(count,1);
                alphaf = cell(count,1);
                model_alphaf = cell(count,1);
                poscount = zeros(count,1);   

            end
            imtut = im;
            issame = 0;
            init = 1;           
        end
         %% of end
		if (size(im,3) > 1) && mod(frame,of_interval) ~= 0,
			im = rgb2gray(im);
		end
		if resize_image && mod(frame,of_interval) ~= 0,
			im = imresize(im, 0.5);
        end
        if isIP == 1,
            if init  == 0 && mod(frame,of_interval) ~= 0,
                issame = psnr(im,imtut) > 40;%max(abs(im(:)-imtut(:))) > 50;%
            elseif mod(frame,of_interval) == 0,
                issame = 0;
            end
        else
            issame = 0;
        end
        if issame == 0 && count ~= 0,
            
            for i = 1:count
                if init  == 0,
                    %obtain a subwindow for detection at the position from last
                    %frame, and convert to Fourier domain (its size is unchanged)
                    patch{i} = get_subwindow(im, pos(i,:), window_sz(i,:));
                    zf{i} = fft2(get_features(patch{i}, features, cell_size,cos_window{i}));%, cos_window{i}
			
                    %calculate response of the classifier at all shifts
                    switch kernel.type
                    case 'gaussian',
                        kzf{i} = gaussian_correlation(zf{i}, model_xf{i}, kernel.sigma);
                    case 'polynomial',
                        kzf{i} = polynomial_correlation(zf{i}, model_xf{i}, kernel.poly_a, kernel.poly_b);
                    case 'linear',
                        kzf{i} = linear_correlation(zf{i}, model_xf{i});
                    end
                    response = real(ifft2(model_alphaf{i} .* kzf{i}));  %equation for fast detection

                    %target location is at the maximum response. we must take into
                    %account the fact that, if the target doesn't move, the peak
                    %will appear at the top-left corner, not at the center (this is
                    %discussed in the paper). the responses wrap around cyclically.
                    peak_value(i) = max(response(:));
                   
                    %% apply new target_sz and pos
                    [vert_delta, horiz_delta] = find(response == peak_value(i), 1);   

                    if vert_delta > size(zf{i},1) / 2,  %wrap around to negative half-space of vertical axis
                        vert_delta = vert_delta - size(zf{i},1);
                    end
                    if horiz_delta > size(zf{i},2) / 2,  %same for horizontal axis
                        horiz_delta = horiz_delta - size(zf{i},2);
                    end                   
                    pos(i,:) = pos(i,:) + cell_size * [vert_delta - 1, horiz_delta - 1];
                end
                
                %obtain a subwindow for training at newly estimated target position
                patch{i} = get_subwindow(im, pos(i,:), window_sz(i,:));
                xf{i} = fft2(get_features(patch{i}, features, cell_size,cos_window{i}));%, cos_window{i}
                
                %Kernel Ridge Regression, calculate alphas (in Fourier domain)
                switch kernel.type
                case 'gaussian',
                    kf{i} = gaussian_correlation(xf{i}, xf{i}, kernel.sigma);
                case 'polynomial',
                    kf{i} = polynomial_correlation(xf{i}, xf{i}, kernel.poly_a, kernel.poly_b);
                case 'linear',
                    kf{i} = linear_correlation(xf{i}, xf{i});
                end
                alphaf{i} = yf{i} ./ (kf{i} + lambda);   %equation for fast training
                %% train
                if init  == 1,  %first frame, train with a single image
                    model_alphaf{i} = alphaf{i};
                    model_xf{i} = xf{i};
                else
                    %subsequent frames, interpolate model
                    model_alphaf{i} = (1 - interp_factor) * model_alphaf{i} + interp_factor * alphaf{i};
                    model_xf{i} = (1 - interp_factor) * model_xf{i} + interp_factor * xf{i};
                end
                %% save position and timing
                positions(i,frame,:) = pos(i,:);
                if init  == 0,
                    if positions(i,frame,:) == positions(i,frame-1,:),
                        poscount(i) = poscount(i) + 1;
                    else
                        poscount(i) = 0;
                    end
                end
            end
       
            
            %% visualisation
            rects = [floor(pos(:,2)) - floor(target_sz(:,2)/2),floor(pos(:,1)) - floor(target_sz(:,1)/2),target_sz(:,2), target_sz(:,1)];

            im_save = insertRect(im,rects,0);

        elseif count == 0;
            im_save = im;
        end

        im_save = drawMask(im_save,roi,0);
        axes(handles.axes1);
        imshow(im_save,'Parent',handles.axes1); 
        drawnow;
        guidata(hObject,handles);
        str = sprintf('%s%d%s',save_path,frame,'.jpg');
        if write_frames || write_video
            imwrite(im_save,str);
        end
        fprintf(fID,'%d %d\n',frame,count);
        fprintf('frame:%d\n',frame);
        frame = frame + 1;
        init = 0;
        imtut = im;
        
    end
    frame = frame - 1;
    set(handles.axes1,'Visible','off');
    set(handles.text_status,'String','Processing the video: Video Finished!');

    fprintf(fID,'fps = %f',time/frame);
    fclose(fID);
    if write_video
        set(handles.text_status,'String','Writing Processed Video!');
        vidWriter = VideoWriter(video_path,'MPEG-4');
        vidWriter.FrameRate = vidReader.FrameRate;
        open(vidWriter);
        for i = 1:frame
            ivideo = imread(sprintf('%s%d%s',save_path,i,'.jpg'));
            writeVideo(vidWriter,ivideo);
        end
        close(vidWriter);
        set(handles.text_status,'String','Processed Video is written, You can replay it!');
        handles.written_video = video_path;
        set(handles.pushbutton_play,'Enable','off'); %% off olacak en son
        set(handles.pushbutton_finish,'Enable','off');
        set(handles.pushbutton_replay,'Enable','on');
        set(handles.edit_video_start,'Enable','off');
        set(handles.edit_video_end,'Enable','off');

    else
        %% clear all data and reset the gui
        clearAll(hObject,handles);
    end

%% end tracker
guidata(hObject,handles);
end


% --- Executes on button press in pushbutton_stop.
function pushbutton_stop_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

% --- Executes on button press in pushbutton_finish.
function pushbutton_finish_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_finish (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global video_finish;
video_finish = 1;

set(handles.pushbutton_finish,'Enable','off');

guidata(hObject,handles);
end

% --- Executes on button press in pushbutton_replay.
function pushbutton_replay_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_replay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global replay_init;
replay_init = ~replay_init;
if replay_init
    video_path = handles.written_video;
    if strcmp(video_path,'') || ~exist(video_path, 'file')
        msgbox('No written video found!','Error','error');
        set(handles.text_status,'String','All data is resetted, please enter from Settings please!');
        set(handles.pushbutton_replay,'String','Finish Playing');
        % clear all data
        return;
    end
    replay_init = ~replay_init;
    vidReader = VideoReader(video_path);
    while hasFrame(vidReader)
        im = readFrame(vidReader);
        set(handles.text_status,'String',sprintf('Playing the Written Video: %.1f seconds',vidReader.CurrentTime));
        axes(handles.axes1);
        imshow(im,'Parent',handles.axes1); 
        drawnow;
    end
else
    msgbox('Written Video has finished playing. All data is resetted, please enter from Settings please!','Warning','warn');
    set(handles.text_status,'String','Written Video has finished playing. All data is resetted, please enter from Settings please!');
    set(handles.pushbutton_replay,'String','Play Written Video');
    %clear all data
    clearAll(hObject,handles);
end
guidata(hObject,handles);
end

function clearAll(hObject,handles)
    %% set gui to beginning state
    cla;
    set(handles.text_status,'String','All data is resetted, please enter from Settings please!');
    set(handles.edit_video_start,'String','');
    set(handles.edit_video_start,'Enable','off');
    set(handles.edit_video_end,'String','');
    set(handles.edit_video_end,'Enable','off');
    set(handles.text_duration,'String','');
    set(handles.checkbox_frames,'Enable','off');
    set(handles.checkbox_video,'Enable','off');
    set(handles.pushbutton_play,'Enable','off');
    set(handles.pushbutton_finish,'Enable','off');
    set(handles.pushbutton_replay,'Enable','off');
    
    %% set handles to beginning states
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
    handles.written_video = '';
    
    %% set globals to beginning states
    global video_finish;
    video_finish = 0;
    global replay_init;
    replay_init = 0;
    
    %% send it to gui
    guidata(hObject,handles);
end

