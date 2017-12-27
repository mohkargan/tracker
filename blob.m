function [count,x,y,width,height] =  blob(FilteredImage)
    [L num]=bwlabel(FilteredImage);
    %min degerler olarak 0.06 secersek ufak arac golgelerini de alabilir,
    %eger 0.07 al?rsak ufak golgeler suzgeclenebilir, eger 0.09 al?rsak tek
    %tek gecen insanlar da sucgeclenebilir
    [image_x,image_y] = size(FilteredImage);

    if (image_x == 540 && image_y == 960) || (image_x == 1080 && image_y) == 1920,
        FRAME_MINY = 0.05;
        FRAME_MINX = 0.05;
    else
        FRAME_MINY = 0.03;
        FRAME_MINX = 0.03;
    end
    FRAME_MAXY = 0.65;
    FRAME_MAXX = 0.65;
    if image_x == 540 && image_y == 960,
        min_area = 400;
    elseif image_x == 1080 && image_y == 1920,
        min_area = 800;
    elseif image_x == 576 && image_y == 1024,
        min_area = 450; 
    elseif (image_x == 288 && image_y == 512) || (image_x == 360 && image_y == 640),
        min_area = 220;   
    elseif image_x == 180 && image_y == 320,
        min_area = 60;   
    elseif image_x == 240 && image_y == 426,
        min_area = 80;
    elseif image_x == 120 && image_y == 213,
        min_area = 40;
    end
    
%     [Bf,Lf,Nf,Af] = bwboundaries(L);
%         
%    figure
%     imshow(L);
% 
%     hold on;
%     for k=1:length(Bf),
% 
%         if(~sum(Af(k,:)))
%             boundary = Bf{k};
%             plot(boundary(:,2), boundary(:,1), 'r','LineWidth',2);
% 
%             for l=find(Af(:,k))'
%                 boundary = Bf{l};
%                 plot(boundary(:,2), boundary(:,1), 'g','LineWidth',2);
%             end
% 
%         end
% 
%     end
    
    STATS=regionprops(L,'all');
    removed=0;
%     Remove the noisy regions
    for i=1:num
        dd=STATS(i).Area;
        if (dd < min_area)
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
%     imshow(L2);
%     %subplot 211,  imshow(L2);%title('BackGround Detected');
%    % subplot 212,  imshow(L2);%title('Blob Detected');
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
    end
end