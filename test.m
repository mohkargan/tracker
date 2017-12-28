close all;clear all;clc;
%% filenames
pathname = 'C:\Users\mohkargan\Desktop\best\video\itü\';
file1 = 'frames_MVI_3159_gaussian_gray_hs_defroi_beginto50.txt';
file2 = 'frames_MVI_3159_gaussian_gray_lk_defroi_beginto50.txt';
file3 = 'frames_MVI_3159_gaussian_hog_hs_defroi_beginto50.txt';
file4 = 'frames_MVI_3159_gaussian_hog_lk_defroi_beginto50.txt';
file5 = 'frames_MVI_3159_polynomial_gray_hs_defroi_beginto50.txt';
file6 = 'frames_MVI_3159_polynomial_gray_lk_defroi_beginto50.txt';
file7 = 'frames_MVI_3159_polynomial_hog_hs_defroi_beginto50.txt';
file8 = 'frames_MVI_3159_polynomial_hog_lk_defroi_beginto50.txt';
file9 = 'frames_MVI_3159_linear_gray_hs_defroi_beginto50.txt';
file10 = 'frames_MVI_3159_linear_gray_lk_defroi_beginto50.txt';
file11 = 'frames_MVI_3159_linear_hog_hs_defroi_beginto50.txt';
file12 = 'frames_MVI_3159_linear_hog_lk_defroi_beginto50.txt';

%% open files to read
fID1 = fopen(fullfile(pathname,file1),'r');
fID2 = fopen(fullfile(pathname,file2),'r');
fID3 = fopen(fullfile(pathname,file3),'r');
fID4 = fopen(fullfile(pathname,file4),'r');
fID5 = fopen(fullfile(pathname,file5),'r');
fID6 = fopen(fullfile(pathname,file6),'r');
fID7 = fopen(fullfile(pathname,file7),'r');
fID8 = fopen(fullfile(pathname,file8),'r');
fID9 = fopen(fullfile(pathname,file9),'r');
fID10 = fopen(fullfile(pathname,file10),'r');
fID11 = fopen(fullfile(pathname,file11),'r');
fID12 = fopen(fullfile(pathname,file12),'r');

%% init values
v1 = zeros(12,298);
v2 = zeros(12,298);
v3 = zeros(12,298);
v4 = zeros(12,298);
v5 = zeros(12,298);
v6 = zeros(12,298);
fID = [fID1, fID2, fID3, fID4, fID5, fID6, fID7, fID8, fID9, fID10, fID11, fID12];

%% assign values
for i = 1:298
    
    for j = 1:12
        tline = fgetl(fID(j)); 
        separated = strsplit(tline,' ');
        v1(j,i) = separated{2}-48;
        v2(j,i) = separated{3}-48;
        v3(j,i) = separated{4}-48;
        v4(j,i) = separated{5}-48;
        v5(j,i) = separated{6}-48;
        v6(j,i) = separated{7}-48;
    end    
  
end

%% ratios
positive_ratio = zeros(12,1);
find_vehicle_ratio = zeros(12,1);
dont_find_ratio = zeros(12,1);
false_find_ratio = zeros(12,1);

for i = 1:12
    positive_ratio(i) = sum(v1(i,:))/sum(v6(i,:));
    find_vehicle_ratio(i) = (sum(v1(i,:))+sum(v4(i,:))+sum(v5(i,:)))/sum(v6(i,:));
    dont_find_ratio(i) = sum(v3(i,:))/sum(v6(i,:));
    false_find_ratio(i) = sum(v2(i,:))/sum(v6(i,:));
end

gauss = zeros(298,1);
lk = zeros(298,1);

for i=1:298
    gauss(i) = (v1(1,i) + v1(3,i) +v1(5,i) + v1(7,i) + v1(9,i) + v1(11,i))/(6*v6(1,i));
    lk(i) = (v1(2,i) + v1(4,i) +v1(1,i) + v1(8,i) + v1(10,i) + v1(12,i))/(6*v6(1,i));
end
plot([1:298], gauss,'r',[1:298], lk,'b')

%% close files
for i = 1:12
    fclose(fID(i));
end