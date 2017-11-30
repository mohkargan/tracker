%  Useful combinations:
%  >> run_tracker choose gaussian hog  %Kernelized Correlation Filter (KCF)
%  >> run_tracker choose linear hog    %Dual Correlation Filter (DCF)
%  >> run_tracker choose gaussian gray %Single-channel KCF (ECCV'12 paper)
%  >> run_tracker choose linear gray   %MOSSE filter (single channel)

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
    
% [filename, pathname] = uigetfile( ...
%    {'*.avi;*.mpg;*.mpeg;*.mp4;*.mkv','Video Files (*.avi,*.mpg,*.mpeg,*.mp4,*.mkv)';
%      '*.*',  'All Files (*.*)'}, ...
%      'Select a video file');

filename = 'MVI_3159.mp4';
pathname = 'C:\Users\mohkargan\Desktop\best\video\itü\';

%call tracker function with all the relevant parameters
[time,framecount] = MTwRCOF(pathname,filename, ...
			padding, kernel, lambda, output_sigma_factor, interp_factor, ...
			cell_size, features);
		
%calculate frames-per-second
fps = framecount / time;
str = sprintf('%12s - FPS:% 4.2f\n', filename,  fps);
        
disp(str)
