%% Test OCR

close all; clc; clear all;
format compact;
fontSize = 20;
addpath("E:\DR_Sonrasi_Projeler\Sempozyumlar\TUFUAB2022\Qgis\cropped_rasters\samples")
tifFiles = dir('E:\DR_Sonrasi_Projeler\Sempozyumlar\TUFUAB2022\Qgis\cropped_rasters\samples\*.tif');
numFiles = length(tifFiles);
mydata = cell(1,numFiles);

%%
% mydata = zeros(numFiles);
for k = 1:numFiles
    mydata{k} = imread(tifFiles(k).name);
end
%%
export = cell(numFiles,2);
%%
for k = 1:numFiles
    % k=1
%     subplot(4,5,k);
    %figure;
    img = mydata{k};

img = imresize(img, 6);

%se = strel('disk',30);
%BW = imtophat(img,se);
%figure
%imshow(BW)
%
BW = rgb2gray(img);
BW = medfilt2(BW, [10 10]);

%imshow(BW)
% 
% Interactively threshold the image.
lowThreshold = 0;
highThreshold = 51;
% highThreshold = 51;

%imhist(BW)
mask = BW < highThreshold;
%imshow(mask)
% Find out areas.
props = regionprops(mask, 'Area');
allAreas = sort([props.Area], 'Descend');

% Keep areas only if they're bigger than 500 pixels and less than 2000 pixels.'
BW = bwareafilt(mask, [500, 20000]);
% n = 1;  
% Idouble = im2double(BW); 
% avg = mean2(Idouble);
% sigma = std2(Idouble);
% BW = imadjust(BW,[avg-n*sigma avg+n*sigma],[]);
% Display the image.
subplot(3, 3, 1);
imshow(img, []);
title('Orijinal Görüntü', 'FontSize', fontSize, 'Interpreter', 'None');
axis('on', 'image');
%hp = impixelinfo();
% Display the image.
subplot(3, 3, 2);
imshow(BW, []);
title('BW Görüntü', 'FontSize', fontSize, 'Interpreter', 'None');
axis('on', 'image');
%hp = impixelinfo();

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
% Get rid of tool bar and pulldown menus that are along top of figure.
% set(gcf, 'Toolbar', 'none', 'Menu', 'none');
% Give a name to the title bar.
set(gcf, 'Name', 'Code: mzeybek', 'NumberTitle', 'on')
drawnow;

%[counts,x] = imhist(BW,16);
%stem(x,counts)
%T = otsuthresh(counts);
%BW = imbinarize(BW,T);


% BW=imbinarize(BW, "adaptive", "Sensitivity", 0.3, "ForegroundPolarity", "dark");
%BW=imbinarize(BW,0.3)

%BW = bwareaopen(BW, 100);

%figure
%bpImage = BW < 150;

%%
% Display the binary image.
subplot(3, 3, 3);
imshow(BW, []);
title('İkili Görüntü', 'FontSize', fontSize, 'Interpreter', 'None');
axis('on', 'image');
%hp = impixelinfo();
drawnow;

%% Remove lines
%  Find branchpoints

bpImage = bwmorph(BW, 'branchpoints'); % Find branch points of skeleton.
bpImage = bwareaopen(bpImage, 20); % Removes all connected components (objects) that have fewer than P pixels from the binary image BW, 
%producing another binary image, bpImage. 
%This operation is known as an area opening.

SE = strel('square', 30); % Morphological structuring element
bpImage2 = imdilate(bpImage, SE); % Dilate image


%%
% Enlarge the branchpoints to really separate the lines.
%se = strel('sphere',6);
%bpImage2 = imdilate(bpImage, 3);
%bpImage2 = imopen(bpImage, ones(1));
linesImage = BW & ~bpImage2;

[labeledImage, numberOfBlobs] = bwlabel(linesImage, 8); % Label connected components in 2-D binary image
% Apply a variety of pseudo-colors to the regions.
coloredLabelsImage = label2rgb (labeledImage, 'hsv', 'k', 'shuffle');
% Display the pseudo-colored image.
subplot(3,3,4)
imshow(coloredLabelsImage);
title('Görüntü Üzerindeki Hatlar', 'FontSize', fontSize, 'Interpreter', 'None');
axis('on', 'image');
%hp = impixelinfo();
drawnow;
%%
% Şimdi her satırdaki pikselleri bulmamız ve bir satıra sığdırmamız gerekiyor
% ve çizginin düz olup olmadığını belirle.
% Çizgi düz değilse, onu tutarız.
% Düz ise, sınır "ızgara" çizgisidir ve biz bunları doldurmak isteyeceğiz.
props = regionprops(linesImage, 'PixelList', 'Solidity'); % Measure properties of image regions
allSolidities = [props.Solidity];
subplot(3,3,5); histogram(allSolidities);
grid on;

linesImage = bwpropfilt(linesImage, 'Solidity', [0.001, 90]); % Extract objects from binary image using properties
% ÖZel noktalarını tekrar yerleştirin.
linesImage = linesImage | bpImage;
% Özel noktaları artık orada olmayan bir çizgiye dokunduysa, bunları kaldırın.
linesImage = bwareaopen(linesImage, 20); % Remove small objects from binary image
% Sekli görüntüleyin.
subplot(3,3,6)
imshow(linesImage);
title('Lines Image', 'FontSize', fontSize, 'Interpreter', 'None');
axis('on', 'image');
hp = impixelinfo();
drawnow;
% Orijinal ikili görüntüyü doldurmak için bu görüntüyü kullanın
%binaryImage2 =~ BW; % Initialize.
%binaryImage2(linesImage) = true;
% Display the image.
%subplot(3,3,7); imshow(binaryImage2);
%title('Final Image', 'FontSize', fontSize, 'Interpreter', 'None');
%axis('on', 'image');
%hp = impixelinfo();
%drawnow;

%skel_image = bwmorph(BW,'skel',3);
%figure
%imshow(skel_image)
%%
% Get rid of blobs touching border
mask = imclearborder(linesImage); % Suppress light structures connected to image border
% Take the 7 largest blobs. 
mask = bwareafilt(mask, 50); % Extract objects from binary image by size
%mask = imrotate(mask, 2);
% Smooth out the numbers by blurring and thresholding.
windowSize = 2;
mask = conv2(double(mask), ones(windowSize)/windowSize^2, 'same') > 0.1; % 2-D convolution
% Display the image.
figure;
imshow(mask, []);
axis('on', 'image');
title('Final Mask', 'FontSize', fontSize, 'Interpreter', 'None');
impixelinfo;
%%
% Rotation
% Detect corners
corners = detectHarrisFeatures(mask); % Detect corners using Harris–Stephens algorithm and return cornerPoints object
 figure;imshow(BW); hold on;
 set(gca, 'ydir', 'reverse')
 plot(corners.selectStrongest(100));
 hold off;
 figure;scatter(corners.Location(:,1),corners.Location(:,2)*-1,"+", "green")

%PCA
pca = pca(corners.Location) % Principal component analysis of raw data
 PcaCol1= sum(pca(1,:)); %the sum of Pca column 1
 PcaCol2= sum(pca(2,:)); %the sum of Pca column 2
%biplot(pca(:,1:2),'Scores',pca(:,1:2)); xlabel("Bileşen 1"); ylabel("Bileşen 2")

bet_a =0;
      NewHigh = max(pca(:)); %getting the highest pca value
%this if-else statement is used to determine the ange of rotation for the image
  if pca(1,1) < 0.9 & pca(2,1) < 0
        %getting the angle of rotation
   bet_a = abs(atan(pca(1,1)/pca(2,1))*180/pi) - 90
     %converts the angle into a negative value
  elseif pca(1,1) > 0.9 & pca(2,1) > 0.2
      bet_a = abs(atan(pca(1,2)/pca(2,2))*180/pi)
  end
  
clear pca
 %BWfilter = bwpropfilt(mask,'perimeter',1);

 %stats = regionprops(BWfilter,'Orientation');  % 56.981 for this image

 % Do the Radon transform.
%theta = 0:90;
%[R,xp] = radon(mask,theta);
% Find the location of the peak of the radon transform image.
%maxR = max(R(:));
%[rowOfMax, columnOfMax] = find(R == maxR)
bet_a
mask = imrotate(BW, bet_a); % Rotate image
%figure;imshow(mask)

% Get the boundaries
[boundaries,LL] = bwboundaries(mask,8,"holes"); % Trace region boundaries in binary image
% Plot them over the inverted gray scale image.
figure ;%imshow(255-grayImage, []);
numBoundaries = length(boundaries);
caption = sprintf('Inverted Grayscale Image with %d boundaries', numBoundaries);
title(caption, 'FontSize', fontSize, 'Interpreter', 'None');
hold on;
for kk = 1 : length(boundaries)
	thisBoundary = boundaries{kk};
	y = thisBoundary(:, 1);
	x = thisBoundary(:, 2);
	plot(x, -y, 'r-', 'LineWidth', 2);
end
%CC = bwconncomp(mask, 4);
Area_S = regionprops(mask, 'Area');
mask=imclearborder(mask,26);
% figure
% imshow(mask)

%% OCR
%ocrResults = ocr(binaryImage2);
ocrResults = ocr(mask); % Recognize text using optical character recognition
%ocrResults = ocr(skel_image);

words = {ocrResults(:).Text}';
words = deblank(words);
words

Iocr= insertObjectAnnotation(img, 'rectangle', ...
                           ocrResults.WordBoundingBoxes, ...
                           ocrResults.WordConfidences);
     figure; imshow(Iocr);
     
%fprintf([k,'.txt'],'%s\n',words{1,1})
fid = fopen("out.txt",'a+');
%fileID = fopen('trytext.txt');
fprintf(fid, '%s\r\n', words{:});
fclose(fid);

end
