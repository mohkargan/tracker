function outIm = drawMask(im,mask,color)
    im_x = size(im,1);
    im_y = size(im,2);
    mask_x = size(mask,1);
    mask_y = size(mask,2);
    outIm = im;
    if (im_x ~= mask_x || im_y ~= mask_y)
        return;
    else       
        if (size(color,2) ~= 3 || size(find(color <= 0 | color>=255),2) ~= 0 )
            color = [255,0,0];
        elseif size(color,2) == 3
            isInt = 1;
            for i=1:3
               if floor(color(i)) ~= color(i)
                   isInt = 0;
                   break;
               end
            end
            if isInt == 0
                color = [255,0,0];
            end
        end
        expanded_mask = zeros(mask_x + 2,mask_y + 2);
        expanded_mask(2:end-1,2:end-1) = mask;
        for i = 2:mask_x + 1
           for j = 2:mask_y + 1   
               left = expanded_mask(i,j-1);
               right = expanded_mask(i,j+1);
               current = expanded_mask(i,j);
               up = expanded_mask(i-1,j);
               down = expanded_mask(i+1,j);
               if size(outIm,3) ~= 3
                   outIm = cat(3, im, im, im);
               end
               if ((left == 0 && current == 1 && right == 1) || (left == 1 && current == 1 && right == 0) || (up == 0 && current == 1 && down == 1) || (up == 1 && current == 1 && down == 0))
                   outIm(i-1,j-1,1) = color(1);
                   outIm(i-1,j-1,2) = color(2);
                   outIm(i-1,j-1,3) = color(3);
               end
           end
        end
    end

end