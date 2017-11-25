function [count,x,y,width,height] =  blob(FilteredImage)
    [L num]=bwlabel(FilteredImage);
    %min degerler olarak 0.06 secersek ufak arac golgelerini de alabilir,
    %eger 0.07 al?rsak ufak golgeler suzgeclenebilir, eger 0.09 al?rsak tek
    %tek gecen insanlar da sucgeclenebilir
    FRAME_MINY = 0.07;
    FRAME_MINX = 0.07;
    FRAME_MAXY = 0.45;
    FRAME_MAXX = 0.45;

    [image_x,image_y] = size(FilteredImage);
    
    STATS=regionprops(L,'all');
    removed=0;
%     Remove the noisy regions
    for i=1:num
        dd=STATS(i).Area;
        if (dd < 800)
            L(L==i)=0;
            removed = removed + 1;
            num=num-1;
        end
    end

    [L2 num2]=bwlabel(L);
    % Trace region boundaries in a binary image.
    %L2 = L;
    [B,L,N,A] = bwboundaries(L2);
    
    %Display results
%     figure
%     subplot 211,  imshow(L2);title('BackGround Detected');
%     subplot 212,  imshow(L2);title('Blob Detected');
% 
%     hold on;
%     for k=1:length(B),
% 
%         if(~sum(A(k,:)))
%             boundary = B{k};
%             plot(boundary(:,2), boundary(:,1), 'r','LineWidth',2);
% 
%             for l=find(A(:,k))'
%                 boundary = B{l};
%                 plot(boundary(:,2), boundary(:,1), 'g','LineWidth',2);
%             end
% 
%         end
% 
%     end
    count = 1;
%     figure 
%     imshow(FilteredImage)
%     figure
%     imshow(L2);
%     hold on;
    for i = 1:length(B)
        boundary = B{i};
        min_x = min(boundary(:,1));
        max_x = max(boundary(:,1));
        min_y = min(boundary(:,2));
        max_y = max(boundary(:,2));
        if ((max_x-min_x)/ image_x < FRAME_MAXX) && ((max_y-min_y)/ image_y < FRAME_MAXY) ...
                && ((max_x-min_x)/ image_x > FRAME_MINX) && ((max_y-min_y)/ image_y > FRAME_MINY)
            x(count) = min_x; %min_x
            width(count) = max_x-min_x; %width
            y(count) = min_y; %min_y
            height(count) = max_y-min_y; %max_y
            
            count = count + 1;
%             rectangle('Position',[min_y min_x (max_y-min_y) (max_x-min_x)],'EdgeColor','g');
        end
    end
%     hold off;
    if count == 1
        x = 0;
        y = 0;
        width = 0;
        height = 0;
%     else
%         x = x(2:end);
%         y = y(2:end);
%         width = width(2:end);
%         height = height(2:end);
    end
end