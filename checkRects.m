function outRects = checkRects(rects)
    count = size(rects,1);
    if count == 1
        outRects = rects;
    else
        roicount = 1;
        inrect = int16.empty;

        new_x = rects(:,1);
        new_y = rects(:,2);
        new_w = rects(:,3);
        new_h = rects(:,4);
        for i = 1:count
            for j = 1:count
                if (i == j)
                    continue;
                end
                if (new_x(j) > new_x(i)) && (new_y(j) > new_y(i)) && (new_x(i)+new_w(i) > new_x(j)+new_w(j))...
                                    && (new_y(i)+new_h(i) > new_y(j)+new_h(j))
                    inrect(roicount) = j;
                    roicount = roicount + 1;
                end
            end
        end
        roicount = roicount - 1;

        if roicount > 0
            outCount = 1;
            outRects = zeros(count-roicount,4);
            for i = 1:count
                I = find(inrect == i);
                hasElement = size(I,2);
                if (hasElement == 0)
                    outRects(outCount,:) = [new_x(i) new_y(i) new_w(i) new_h(i)];
                    outCount = outCount +1;
                end            
            end
        else 
            outRects = rects;
        end
    end
end