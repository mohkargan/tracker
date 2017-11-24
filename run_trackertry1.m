%  Useful combinations:
%  >> run_tracker choose gaussian hog  %Kernelized Correlation Filter (KCF)
%  >> run_tracker choose linear hog    %Dual Correlation Filter (DCF)
%  >> run_tracker choose gaussian gray %Single-channel KCF (ECCV'12 paper)
%  >> run_tracker choose linear gray   %MOSSE filter (single channel)

%function [precision, fps] = run_tracker(video, kernel_type, feature_type, show_visualization, show_plots)
close all;  clear all; clc;

    kernel_type = 'gaussian';
    feature_type = 'hog';
    
	kernel.type = kernel_type;
	
	features.gray = false;
	features.hog = false;
	
	padding = 1.5;  %extra area surrounding the target
	lambda = 1e-4;  %regularization
	output_sigma_factor = 0.1;  %spatial bandwidth (proportional to target)
	
	switch feature_type
	case 'gray',
		interp_factor = 0.075;  %linear interpolation factor for adaptation

		kernel.sigma = 0.2;  %gaussian kernel bandwidth
		
		kernel.poly_a = 1;  %polynomial kernel additive term
		kernel.poly_b = 7;  %polynomial kernel exponent
	
		features.gray = true;
		cell_size = 1;
		
	case 'hog',
		interp_factor = 0.02;
		
		kernel.sigma = 0.5;
		
		kernel.poly_a = 1;
		kernel.poly_b = 9;
		
		features.hog = true;
		features.hog_orientations = 9;
		cell_size = 4;
		
	otherwise
		error('Unknown feature.')
	end


	assert(any(strcmp(kernel_type, {'linear', 'polynomial', 'gaussian'})), 'Unknown kernel.')
    
    [filename, pathname] = uigetfile( ...
    {'*.avi;*.mpg;*.mpeg;*.mp4;*.mkv','Video Files (*.avi,*.mpg,*.mpeg,*.mp4,*.mkv)';
     '*.*',  'All Files (*.*)'}, ...
     'Select a video file');

    fullfilename = fullfile(pathname,filename);

    vidReader = VideoReader(fullfilename);

    framecount = 1;
    init = 1;
    time = tic();
    while hasFrame(vidReader)
        im = readFrame(vidReader);
        if framecount == 2
            init = 0;
        end
        if framecount == 1
            imshow(im);
            rect = getrect;
            pos = [floor(rect(1)) floor(rect(2))];
            target_sz = [floor(rect(3)) floor(rect(4))];
            model_alphaf = 0;
            model_xf = 0;
            
            resize_image = (sqrt(prod(target_sz)) >= 100);  %diagonal size >= threshold
            if resize_image,
                pos = floor(pos / 2);
                target_sz = floor(target_sz / 2);
            end


            %window size, taking padding into account
            window_sz = floor(target_sz * (1 + padding));

	
            %create regression labels, gaussian shaped, with a bandwidth
            %proportional to target size
            output_sigma = sqrt(prod(target_sz)) * output_sigma_factor / cell_size;
            yf = fft2(gaussian_shaped_labels(output_sigma, floor(window_sz / cell_size)));

            %store pre-computed cosine window
            cos_window = hann(size(yf,1)) * hann(size(yf,2))';	          
            
        end
        %call tracker function with all the relevant parameters
        [pos2,target_sz2,model_alphaf2,model_xf2] = trackertry1(im,init,model_alphaf,model_xf,resize_image,window_sz,yf,cos_window, pos,...
            target_sz, kernel, lambda, interp_factor, cell_size, features);
        pos = pos2;
        target_sz = target_sz2;
        model_alphaf = model_alphaf2;
        model_xf = model_xf2;
        framecount = framecount + 1;
        imshow(imresize(im, 0.5));
        rectangle('Position',[pos , target_sz([2,1])],'EdgeColor','g');%- target_sz([2,1])/2
        %pos = floor(pos - target_sz([2,1])/2);
        disp('pos ');disp(pos);
        disp('target_sz ');disp(target_sz);
    end
    time = time + toc();
	%calculate frames-per-second
	fps = framecount / time;

    str = sprintf('%12s - FPS:% 4.2f\n', filename,  fps);
        
    disp(str)
