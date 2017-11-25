function [count,x,y,width,height] = ofmod(opticFlow, frameGray1,frameGray4)
    flow = estimateFlow(opticFlow,frameGray1); 
%     flow = estimateFlow(opticFlow,frameGray2); 
%     flow = estimateFlow(opticFlow,frameGray3); 
    flow = estimateFlow(opticFlow,frameGray4); 

    fm = flow.Magnitude;
    bfm = fm > 0.01;
    frameMedian = medfilt2(bfm);        % median filter
    se = strel('disk',45);%high resolution için 45, dü?ük için 15
    closeBW = imclose(frameMedian,se);  % morphological close
    se1 = strel('line',11,90);
    erodedBW = imerode(closeBW,se1);    % erode
    [count,x,y,width,height] = blob(erodedBW);
    count = count - 1;
end