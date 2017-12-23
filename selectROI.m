function mask = selectROI(I)
    BW = roipoly(I);
    mask = BW;
    imshow(mask);
end