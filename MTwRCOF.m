% MultiTracker with Regularly Checked Optical Flow
function [time,frame] = MTwRCOF(pathname,filename,  ...
	padding, kernel, lambda, output_sigma_factor, interp_factor,...
    cell_size, features, defaultroi, start_time, isIP, opticalFlow)

    %% open video file
    fullfilename = fullfile(pathname,filename);
    vidReader = VideoReader(fullfilename);
    
    %% find total frame count
    framecount = floor(vidReader.Duration * vidReader.FrameRate) + 1;
    
    %% set starting time
    if start_time >= 0 && vidReader.Duration
        vidReader.CurrentTime = start_time;
    end

    %% set roi
    if defaultroi == 1
        roistruct = load('MVI_mask.mat','-mat');
        roi = roistruct.BW;
    else
        %kullanýcýnýn alaný sol ile seçip, birleþtirip en son alana sað
        %týklatýp "create mask" demesi gerekmektedir
        im = readFrame(vidReader);
        im = imresize(im,0.5);
        roi = roipoly(im);
    end
    
    %% set saving path according to folder of video, correlation, features and
    filename_noext = strtok(filename,'.'); %filename without extention
    save_path = sprintf('%s%s%s',pathname,'frames_',filename_noext);
    
    save_path = strcat(save_path,'_',kernel.type,'_');
    
    if features.hog 
        save_path = strcat(save_path,'hog');
    else
        save_path = strcat(save_path,'gray');
    end
    
    if defaultroi == 1
        save_path = strcat(save_path,'_defroi\'); %default roi
    else
        save_path = strcat(save_path,'_userroi\'); % user selected roi
    end
        
    mkdir(save_path);

    %% optical flow ile otomatik seçme   
    if strcmp(opticalFlow.type,'Lucas-Kanade')
        opticFlow = opticalFlowLK('NoiseThreshold',opticalFlow.noiseThreshold);
    else %if strcmp(opticalFlow.type,'Horn-Schunck')
        opticFlow = opticalFlowHS('MaxIteration',10,'Smoothness',1);
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
	time = 0;  %to calculate FPS
    
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

	while hasFrame(vidReader),
%         if frame == 75
%             stop = 1;
%         end
        %% 20 frame boyunca ayni pozisyondaysa trackerdan sil
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
        if count ~=0,
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
                opticFlow = opticalFlowLK('NoiseThreshold',opticalFlow.noise_threshold);
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
            tic()
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
        
            time = time + toc();
            
            %% visualisation
            rects = [floor(pos(:,2)) - floor(target_sz(:,2)/2),floor(pos(:,1)) - floor(target_sz(:,1)/2),target_sz(:,2), target_sz(:,1)];

            im_save = insertRect(im,rects,0);

        elseif count == 0;
            im_save = im;
        end

        im_save = drawMask(im_save,roi,0);
        imshow(im_save)
        str = sprintf('%s%d%s',save_path,frame,'.jpg');
        imwrite(im_save,str);
            
        fprintf('frame:%d\n',frame);
        frame = frame + 1;
        init = 0;
        imtut = im;

    end
    frame = frame - 1;
end

