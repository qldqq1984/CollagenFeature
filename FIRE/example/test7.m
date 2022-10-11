clear all;
close all;

tic
addpath(genpath('../')); 
addpath('e:\otsu');
mfn = mfilename;
%Img = imread('d:\NSP.tif');
%green = Img(:,:,2);

%height = size(green,1);
%width = size(green,2);

%load('d:\Anju\NCIC-1.mat');
raw_pic_name = 'e:\otsu\tilescan.lsm';
raw = tiffread(raw_pic_name);
%raw = imread(raw_pic_name);
green = raw.green;
height = raw.height;
width = raw.width;
%green = raw(:,:,2);
%height = size(raw,1);
%width = size(raw,2);

%[mask,a,b,c] = EMSeg(green,2);
%mask = mask - 1;
mask = uint16(otsu(green,2));
se = strel('diamond',1);
mask = imerode(mask,se);
mask = imdilate(mask,se);
%imshow(mask);


p.Nimages   = 1;  
p.yred      = 1:height;
p.xred      = 1:width; 
p = param_example(p); 


im3 = zeros(1,height,width);
im3(1,:,:) = green.*uint16(mask);
data = fire(p,im3,2);
%figure(4);
%imagesc(mask);

%disp(degree1);


fiber = data.F;
points = data.X;
fiber_degree = zeros(size(fiber,2),1);
fiber_length = zeros(size(fiber,2),1);
%diff_degree = zeros(size(fiber,2),1);

for k=1:1:size(fiber,2)
    temp = fiber(k).v;
    tx = zeros(size(temp,2),1);
    ty = zeros(size(temp,2),1);
    temp_length = 0;
    for j=1:1:size(temp,2)
        tx(j) = points(temp(j),1);
        ty(j) = points(temp(j),2);
        if j~=1
            temp_length = temp_length + sqrt( (tx(j)-tx(j-1))^2 + (ty(j)-ty(j-1))^2 );
        end
    end
    temp_p = polyfit(tx,ty,1);
    fiber_degree(k,1) = atand(temp_p(1));
    fiber_length(k,1) = temp_length;
    %diff_degree(k,1) = abs(degree1 - fiber_degree(k,1));
    %if diff_degree(k,1) > 90
    %    diff_degree(k,1) = 180-diff_degree(k,1);
    %end
end


figure(3);
hist(fiber_degree);
figure(4);
hist(fiber_length);
%}


