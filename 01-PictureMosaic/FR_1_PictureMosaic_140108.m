%===================================================
%
% Picture Mosaic
%
% Author: Chayut Orapinpatipat 2013
% 
%===================================================

%%
clear all;
close all;
clc;
%%
%---------User Defined Varibles-------%

nLibsize = 77;   %number of images in library
nCellrow = 64;   
nCellcolumn = 64;
Cellwidth = 81;
Cellheight = 115;
FAK = 0.62; %Weighing factor of small images
OS =10;     %Image Brightness offset
%C:\\Users\\Chayut\\Dropbox\\Academic File 3.2\\Intern (img)
FNAMEFMT = '..\\..\\..\\Resource\\project 1 numbered img\\1 (%d).jpg'; %Images Library
FNAMEOUT2 = 'Kleinmann_FAK%d_ScrHisteq_OS%d.jpg'; %Name for output file

A =imread('Kleinmann.jpg'); %Template Image File

%%
Baseheight = nCellrow * Cellheight;
Basewidth = nCellcolumn * Cellwidth;
IMAGES = cell(1,nLibsize);
GRAYIMAGES = cell(1,nLibsize);
BWaverage = zeros(nLibsize);
FNAMEOUT2 = sprintf(FNAMEOUT2, FAK*100, OS);

%%
% Load images
for i=1:nLibsize
  IMAGES{i} = imread(sprintf(FNAMEFMT, i)); %load image into cell
  GRAYIMAGES{i}= rgb2gray(IMAGES{i}) ;     %convert RGB image to Gray
  BWaverage(i) = uint8(mean2(IMAGES{i})); %find average of each picture
end
%================================================

Agray = rgb2gray(A);
Abig = imresize(Agray,[Baseheight,Basewidth]); %resize input template image to final image size

B = uint8(zeros(Baseheight,Basewidth)); %allocaltion for output image

%loop for each cell
for i = 0:1:(nCellrow-1)
    for j = 0:1:(nCellcolumn-1)  
        n =randi([1, nLibsize]); %randomized image selection for each cell 
        source1 = GRAYIMAGES{n}; %retrive image for library
        average1 = BWaverage(n); %obtain the image averaged brightness
        %loop for each pixel in each pixel
        for k = 1:1:Cellwidth
            for m = 1:1:Cellheight
                temp = (double(Abig(m+(Cellheight*i),k+(Cellwidth*j)))-127)*FAK +127 + (double(source1(m,k)) - average1) * (1-FAK) +OS; %combination of pixel values
                B(m+(Cellheight*i),k+(Cellwidth*j))=   uint8(temp); %convert to uint8 then place in desinated location
            end
        end
    end
end


imwrite(B,FNAMEOUT2)
figure(1)
imshow(B);
figure(2)
imhist(B);
 


