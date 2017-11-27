function img_ret = insertRect(grayImage,rects,color)
    img = cat(3, grayImage, grayImage, grayImage);
    [img_x img_y] = size(grayImage);
    if color == 0,
%         color = randi([0 25],1,3);
        color = [0 255 0];
    end
    for i = 1:size(rects,1),
       x = rects(i,2); y = rects(i,1); w = rects(i,4); h = rects(i,3);
       if (x+w < img_x) && (w < img_x)&& (x < img_x)&& (y < img_y) && (y+h < img_y) && (h < img_y) && (x > 1) && (y > 1) && (w > 1) && (h > 1),
        for j= 1:3,
            img(x:x+w,y,j)   = color(j);
            img(x:x+w,y+h,j) = color(j);
            img(x,y:y+h,j)   = color(j);
            img(x+w,y:y+h,j) = color(j);   
       end
       end
    end
    img_ret = img;
end