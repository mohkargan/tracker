function [count,x,y,width,height] = ofmod(opticFlow, frameGray)

    flow = estimateFlow(opticFlow,frameGray); 

    fm = flow.Magnitude;
    bfm = fm > 0.01;
    frameMedian = medfilt2(bfm);        % median filter
    se = strel('disk',15);
    closeBW = imclose(frameMedian,se);  % morphological close
    se1 = strel('line',11,90);
    erodedBW = imerode(closeBW,se1);    % erode
    [count,x,y,width,height] = blob(erodedBW);
    count = count - 1;
end