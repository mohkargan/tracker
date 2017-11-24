%function [positions, time] = tracker(video_path, img_files, pos, target_sz, ...
%	padding, kernel, lambda, output_sigma_factor, interp_factor, cell_size, ...
%	features, show_visualization)
function [pos,target_sz, model_alphaf, model_xf] = trackertry1(im,init,model_alphaf,model_xf,resize_image,window_sz,...
    yf,cos_window,pos,...
    target_sz, kernel, lambda, interp_factor, cell_size,features)
	%if the target is large, lower the resolution, we don't need that much
	%detail
% 	resize_image = (sqrt(prod(target_sz)) >= 100);  %diagonal size >= threshold
% 	if resize_image,
% 		pos = floor(pos / 2);
% 		target_sz = floor(target_sz / 2);
% 	end
% 
% 
% 	%window size, taking padding into account
% 	window_sz = floor(target_sz * (1 + padding));
% 	
% % 	%we could choose a size that is a power of two, for better FFT
% % 	%performance. in practice it is slower, due to the larger window size.
% % 	window_sz = 2 .^ nextpow2(window_sz);
% 
% 	
% 	%create regression labels, gaussian shaped, with a bandwidth
% 	%proportional to target size
% 	output_sigma = sqrt(prod(target_sz)) * output_sigma_factor / cell_size;
% 	yf = fft2(gaussian_shaped_labels(output_sigma, floor(window_sz / cell_size)));
% 
% 	%store pre-computed cosine window
% 	cos_window = hann(size(yf,1)) * hann(size(yf,2))';	
	
	%note: variables ending with 'f' are in the Fourier domain.

	if size(im,3) > 1,
		im = rgb2gray(im);
	end
	if resize_image,
		im = imresize(im, 0.5);
	end

	if ~init,
		%obtain a subwindow for detection at the position from last
		%frame, and convert to Fourier domain (its size is unchanged)
		patch = get_subwindow(im, pos, window_sz);
		zf = fft2(get_features(patch, features, cell_size, cos_window));
		
		%calculate response of the classifier at all shifts
		switch kernel.type
		case 'gaussian',
			kzf = gaussian_correlation(zf, model_xf, kernel.sigma);
		case 'polynomial',
			kzf = polynomial_correlation(zf, model_xf, kernel.poly_a, kernel.poly_b);
		case 'linear',
			kzf = linear_correlation(zf, model_xf);
		end
		response = real(ifft2(model_alphaf .* kzf));  %equation for fast detection

		%target location is at the maximum response. we must take into
		%account the fact that, if the target doesn't move, the peak
		%will appear at the top-left corner, not at the center (this is
		%discussed in the paper). the responses wrap around cyclically.
		[vert_delta, horiz_delta] = find(response == max(response(:)), 1);
		if vert_delta > size(zf,1) / 2,  %wrap around to negative half-space of vertical axis
			vert_delta = vert_delta - size(zf,1);
		end
		if horiz_delta > size(zf,2) / 2,  %same for horizontal axis
			horiz_delta = horiz_delta - size(zf,2);
		end
		pos = pos + cell_size * [vert_delta - 1, horiz_delta - 1];
	end

	%obtain a subwindow for training at newly estimated target position
	patch = get_subwindow(im, pos, window_sz);
	xf = fft2(get_features(patch, features, cell_size, cos_window));

	%Kernel Ridge Regression, calculate alphas (in Fourier domain)
	switch kernel.type
	case 'gaussian',
		kf = gaussian_correlation(xf, xf, kernel.sigma);
	case 'polynomial',
		kf = polynomial_correlation(xf, xf, kernel.poly_a, kernel.poly_b);
	case 'linear',
		kf = linear_correlation(xf, xf);
	end
	alphaf = yf ./ (kf + lambda);   %equation for fast training

    if init  %first frame, train with a single image
		model_alphaf = alphaf;
		model_xf = xf;
    else
		%subsequent frames, interpolate model
		model_alphaf = (1 - interp_factor) * model_alphaf + interp_factor * alphaf;
		model_xf = (1 - interp_factor) * model_xf + interp_factor * xf;
    end
%     box = [pos([2,1]) - target_sz([2,1])/2, target_sz([2,1])];%pos([2,1]) - target_sz([2,1])/2
%     box = [floor(box(1)) floor(box(2)) floor(box(3)) floor(box(4))];
%     disp(box)
end