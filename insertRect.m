function img_ret = insertRect(grayImage,rects,color)
    if size(grayImage,3) ~= 3
        img = cat(3, grayImage, grayImage, grayImage);
    end
    [img_x, img_y] = size(grayImage);
    if (size(color,2) ~= 3 || size(find(color <= 0 | color>=255),2) ~= 0 )
    	color = [0,255,0];
    elseif size(color,2) == 3
    	isInt = 1;
        for i=1:3
        	if (floor(color(i)) ~= color(i))
            	isInt = 0;
                break;
            end
        end
        if isInt == 0
        	color = [255,0,0];
        end
    end
    for i = 1:size(rects,1),
       x = rects(i,2); y = rects(i,1); w = rects(i,4); h = rects(i,3);
       if x+w == img_x
           w = w - 1;
       end
       if y + h == img_y
           h = h - 1;
       end
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