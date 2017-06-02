clear all;
close all;
clc;

% ************************************************************************
% Kris Cote
% BEE 515
% Final Project: Power line removal using matched filters
% ************************************************************************

img=imread('powerlines.jpg');

im=img(:,:,2);                            	% select green plane

figure(1)                                   % show original image and green plane
subplot(1,2,1)
imshow(img,[])
title('Original image')
subplot(1,2,2)
imshow(im,[])
title('Green plane')


%  *************************  Segmentation  *******************************

kernArray=load('centeredKern.mat');
kernArray=double(kernArray.centeredKern);

%normalize kern arrays
for i=1:size(kernArray,3)
    kernArray(:,:,i)=kernArray(:,:,i).*(10/max(max(kernArray(:,:,i))));    
end

kernArray=double(kernArray);

figure(2)
imshow(kernArray(:,:,5))
title('Kernel')

%  **********    Initial Filter
h=fspecial('log',[10,10],1);
imBlurred=imfilter(im,h);


%  **********    Matched Filters
n=12;
temp=abs(conv2(imBlurred,kernArray(:,:,1)));
mfr=zeros(size(temp,1),size(temp,2),n);
mfr(:,:,1)=temp/max(max(temp));
avgVec=zeros(n,1);
avgVec(1)=sum(sum(mfr(:,:,1)))/(size(mfr,1)*size(mfr,2));

for i=2:n
    temp=abs(conv2(imBlurred,kernArray(:,:,i)));
    mfr(:,:,i)=temp;
    avgVec(i)=mean(mean(temp));
end

% normalize across mfr
minAve=max(avgVec);
for i=1:n
    temp=mfr(:,:,i)*(minAve/avgVec(i));
    mfr(:,:,i)=temp/max(max(temp));
end

%select maximum value for each indicies
mfrResult=zeros(size(mfr,1),size(mfr,2));
for i=1:size(mfr,1)
    for j=1:size(mfr,2)
        mfrResult(i,j)=max(mfr(i,j,:));
    end
end



% **********  NORMALIZE MfrRESULTS

resMin=min(min(mfrResult));
mfrResult=mfrResult+abs(resMin);
resMax=max(max(mfrResult));
mfrResult=1-mfrResult/resMax;

mfrResult=1-mfrResult;                  % invert result of convolution
                                        % to make background dark

figure(3)
subplot(1,2,1)
imshow(mfr(:,:,5),[])
title('One matched filter result')
subplot(1,2,2)
imshow(mfrResult,[])
title('Combined matched filter result')

% *********** GET HISTOGRAM

imMax=max(max(mfrResult));
Nbins=255;
hist=zeros(Nbins+1,1);
for i=1:size(mfrResult,1)
    for j=1:size(mfrResult,2)   
        index=ceil((mfrResult(i,j)/imMax)*Nbins)+1;
        hist(index)=hist(index)+1;
    end
end

figure(4)
plot(hist)
title('Histogram of Matched Filter Result')
xlabel('Pixel Value')
ylabel('Pixel Count')

% *********** Find thresholding value starting at maximum bin count

[binMax,Maxi]=max(hist);
aveVal=1;
minAve=0.015;
thold = Maxi/255;
while aveVal>minAve               
    tempImage=im2bw(mfrResult,thold);
    aveVal=sum(sum(tempImage))/(size(tempImage,1)*size(tempImage,2));
    if aveVal <= minAve
       thresholdImage = tempImage;
    end
    thold=thold+0.001;
end

% ******* Morphology on binary threshold image

% Dilation
se=strel('disk',2,0);
dilated1=imdilate(thresholdImage,se);

% Skeleton
skeleton=uint8(bwmorph(dilated1,'skel',Inf));

% Dilation 
se=strel('disk',5,0);
skelSize=size(skeleton);
dilated2=uint8(imdilate(skeleton,se));

figure(5)
subplot(1,3,1)
imshow(dilated1,[])
title('Dilated')

subplot(1,3,2)
imshow(skeleton,[])
title('Skeleton')

subplot(1,3,3)
imshow(dilated2,[])
title('Re-Dilated')

segmentMask=dilated2;

%  *************************  INTERPOLATOIN *******************************

imSize=size(img);
output=uint8(zeros(imSize));
for plane=1:3
    im=img(:,:,plane);
    im=im.*(254/255);                   % scale image down by 1 intensity value
                                        % so that we can differentiate the
                                        % background from the dilated threshold
                                        % mask (valued at 255)
    overlay=uint8(im)+255.*segmentMask((skelSize(1)-imSize(1))/2:skelSize(1)+(imSize(1)-skelSize(1))/2-1,...
        (skelSize(2)-imSize(2))/2:imSize(2)+(skelSize(2)-imSize(2))/2-1);



    region=zeros(5,5);
    pad=[2,2];

    % pad overlay array to avoid edge effects
    overlay=padarray(overlay,pad,'both');

    % guess the contents of masked pixels based on content of neighbors
    for i=1:imSize(1)
        for j=1:imSize(2)
            if overlay(i,j)==255
                pxlCnt=0.0;
                pxlSum=0;
                pxlSum=double(pxlSum);
                % interpolation: average of neighbors with value <255
                region=overlay(i-2:i+2,j-2:j+2);
                    for x=1:5
                        for y=1:5
                            if region(x,y) < 255
                                pxlCnt = pxlCnt+1;
                                pxlSum = pxlSum+double(region(x,y));
                            end
                        end
                    end
                overlay(i,j)=uint8(pxlSum/pxlCnt);
            end
        end
    end
    output(:,:,plane)=overlay(2:imSize(1)+1,2:imSize(2)+1);
end

figure(6)
imshow(output,[])
title('Final image')




