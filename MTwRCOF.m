% MultiTracker with Regularly Checked Optical Flow
function [time,frame] = MTwRCOF(pathname,filename,  ...
	padding, kernel, lambda, output_sigma_factor, interp_factor,...
    cell_size, features)

    fullfilename = fullfile(pathname,filename);

    vidReader = VideoReader(fullfilename);
    framecount = vidReader.NumberOfFrames;
    vidReader = VideoReader(fullfilename);
    
    save_path = sprintf('%s%s%s/',pathname,'/frames_',filename);
    mkdir(save_path);
    isIP = 0;
    %% optical flow ile otomatik seçme   
        %% with ofmod func
    opticFlow = opticalFlowHS;
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

    [count,x,y,width,height] = ofmod(opticFlow,frameGray1, frameGray);
%     cx = x+width/2;
%     cy = y+height/2;
%     %% ROI 
%     roicount = 1;
%     for i = 1:count,
%         if (cy(i) > 115) && (cy(i) >752.5 - 1.77*cx(i)) && (cy(i) <1.44*cx(i)-893),
%         	inroi(roicount) = i;
%             roicount = roicount + 1;
%         end
%     end
%     count = roicount - 1;
%     if count ~= 0,
%      	target_sz = [width(inroi), height(inroi)];
%         pos = [cx(inroi) cy(inroi)];
%     else
%         target_sz = int16.empty;
%         pos = int16.empty;
%     end
    %% merkezi ve pencere buyuklugunu belirle
    target_sz= zeros(count,2);
    pos = zeros(count,2);
    for i=1:count
        target_sz(i,:) = [width(i), height(i)];
        pos(i,:) = [x(i), y(i)]+ floor(target_sz(i,:)/2);%merkez
    end
    
    resize_image = 1;
%     for i=1:count
%         if (sqrt(prod(target_sz)) >= 100) == 1
%             resize_image = 1;
%             break;
%         end
%     end
    
	if resize_image && ~isempty(pos),
		pos = floor(pos / 2);
		target_sz = floor(target_sz / 2);
    end
    
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
    get_small = 0;
    get_big = 0;
    frame = 1;
    issame = 0;
    init = 1;
    of_interval = 15;
    del_interval = 10;
    scale_step = 1.05;
    scale_weight = 0.95;  
    

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
        %% of buraya
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
            opticFlow = opticalFlowHS;
            [count,x,y,width,height] = ofmod(opticFlow,imtut, im);
            cx = x+width/2;
            cy = y+height/2;
            %% ROI 
            roicount = 1;
            for i = 1:count,
                if (cy(i) > 115) && (cx(i) >752.5 - 1.77*cy(i)) && (cx(i) <1.44*cy(i)-663),
                    inroi(roicount) = i;
                    roicount = roicount + 1;
                end
            end
            count = roicount - 1;
            if count ~= 0,
            target_sz = [width(inroi), height(inroi)];
            pos = [cx(inroi) cy(inroi)];
            
            %% merkezi ve pencere buyuklugunu belirle
%             target_sz= zeros(count,2);
%             pos = zeros(count,2);
%             for i=1:count
%                 target_sz(i,:) = [width(i), height(i)];
%                 pos(i,:) = [x(i), y(i)]+ floor(target_sz(i,:)/2);%merkez
%             end
    
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
                    
                    
                    %% target_sz adjustment
%                     new_window_sz = floor(window_sz(i,:)./scale_step);
%                     
%                     output_sigma = sqrt(prod(floor(target_sz(i,:).*scale_step))) * output_sigma_factor / cell_size;
%                     new_yf = fft2(gaussian_shaped_labels(output_sigma, floor(new_window_sz / cell_size)));
% 
%                     new_cos_window = hann(size(new_yf,1)) * hann(size(new_yf,2))';	
%                     
%                     new_patch = get_subwindow(im, pos(i,:), new_window_sz);
%                     new_zf = fft2(get_features(new_patch, features, cell_size,new_cos_window));
%                     new_model_xf = zeros(size(new_zf));
%                     for j = 1:size(new_zf,3)
%                         new_model_xf(:,:,j) = imresize(model_xf{i}(:,:,j),[size(new_zf,1) size(new_zf,2)]);
%                     end
%                     
%                     switch kernel.type
%                     case 'gaussian',
%                         new_kzf = gaussian_correlation(new_zf, new_model_xf, kernel.sigma);
%                     case 'polynomial',
%                         new_kzf = polynomial_correlation(new_zf, new_model_xf, kernel.poly_a, kernel.poly_b);
%                     case 'linear',
%                         new_kzf = linear_correlation(new_zf, new_model_xf);
%                     end
%                                         
%                     new_model_alphaf = imresize(model_alphaf{i},[size(new_zf,1) size(new_zf,2)]);
%                     
%                     new_res = real(ifft2(new_model_alphaf .* new_kzf));
%                     new_peak_value = max(new_res(:));
%                     
%                     if (scale_weight * new_peak_value > peak_value(i)),
%                         target_sz(i,:) = floor(target_sz(i,:)/scale_step);
%                         get_big = get_big + 1;
%                         response = new_res;
%                         peak_value(i) = new_peak_value;
%                         
%                         window_sz(i,:) = new_window_sz;
%                         yf{i} = new_yf;
%                         cos_window{i} = new_cos_window;
%                         xf{i} = new_xf;
%                         kf{i} = new_kf;
%                         model_alphaf{i} = new_model_alphaf;
%                         model_xf{i} = new_model_xf;                       
%                     end
%                     
%                     new_window_sz = floor(window_sz(i,:).*scale_step);
%                     
%                     output_sigma = sqrt(prod(floor(target_sz(i,:).*scale_step))) * output_sigma_factor / cell_size;
%                     new_yf = fft2(gaussian_shaped_labels(output_sigma, floor(new_window_sz / cell_size)));
% 
%                     new_cos_window = hann(size(new_yf,1)) * hann(size(new_yf,2))';	
%                     
%                     new_patch = get_subwindow(im, pos(i,:), new_window_sz);
%                     new_zf = fft2(get_features(new_patch, features, cell_size,new_cos_window));
%                     new_model_xf = zeros(size(new_zf));
%                     for j = 1:size(new_zf,3)
%                         new_model_xf(:,:,j) = imresize(model_xf{i}(:,:,j),[size(new_zf,1) size(new_zf,2)]);
%                     end
%                     
%                     switch kernel.type
%                     case 'gaussian',
%                         new_kzf = gaussian_correlation(new_zf, new_model_xf, kernel.sigma);
%                     case 'polynomial',
%                         new_kzf = polynomial_correlation(new_zf, new_model_xf, kernel.poly_a, kernel.poly_b);
%                     case 'linear',
%                         new_kzf = linear_correlation(new_zf, new_model_xf);
%                     end
%                                         
%                     new_model_alphaf = imresize(model_alphaf{i},[size(new_zf,1) size(new_zf,2)]);
%                     
%                     new_res = real(ifft2(new_model_alphaf .* new_kzf));
%                     new_peak_value = max(new_res(:));
% 
%                     if (scale_weight * new_peak_value > peak_value(i)),
%                         target_sz(i,:) = floor(target_sz(i,:)*scale_step);
%                         get_small = get_small + 1;
%                         response = new_res;
%                         peak_value(i) = new_peak_value;
%                         
%                         window_sz(i,:) = new_window_sz;
%                         yf{i} = new_yf;
%                         cos_window{i} = new_cos_window;
%                         xf{i} = new_xf;
%                         kf{i} = new_kf;
%                         model_alphaf{i} = new_model_alphaf;
%                         model_xf{i} = new_model_xf;  
%                     end
                    
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

