function [count,x,y,width,height] = ofmod(opticFlow, frameGray1,frameGray2)
    flow = estimateFlow(opticFlow,frameGray1);  
    flow = estimateFlow(opticFlow,frameGray2); 
    [image_x,image_y] = size(frameGray1);

    fm = flow.Magnitude;
    bfm = fm > 0.01;
    frameMedian = medfilt2(bfm);        % median filter
    if (image_x == 540 && image_y == 960),
        disk_r = 15;
    elseif image_x == 1080 && image_y == 1920,
        disk_r = 45;
    elseif image_x == 576 && image_y == 1024,
        disk_r = 20;
    elseif (image_x == 288 && image_y == 512) || (image_x == 360 && image_y == 640),
        disk_r = 12;
    elseif image_x == 240 && image_y == 426,
        disk_r = 9;
    elseif image_x == 180 && image_y == 320,
        disk_r = 7;
    elseif image_x == 120 && image_y == 213,
        disk_r = 5;
    end
    se = strel('disk',disk_r);
    closeBW = imclose(frameMedian,se);  % morphological close
    se1 = strel('line',11,90);
    erodedBW = imerode(closeBW,se1);    % erode
            
%         figure
%         imshow(frameGray2)
        
%         figure
%         imshow(frameGray2) 
%         hold on
%         plot(flow,'DecimationFactor',[5 5],'ScaleFactor',20)
%         hold off
%         
%         figure
%         imshow(bfm)
%         
%         figure
%         imshow(frameMedian)
%         
%         figure
%         imshow(closeBW)
%         
%         figure
%         imshow(erodedBW)
    
    [count,x,y,width,height] = blob(erodedBW);
    count = count - 1;
end