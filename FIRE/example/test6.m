clear all;
close all;

tic
addpath(genpath('../')); 
mfn = mfilename;
%Img = imread('d:\NSP.tif');
%green = Img(:,:,2);

%height = size(green,1);
%width = size(green,2);

%load('d:\Anju\NCIC-1.mat');
raw_pic_name = 'd:\Anju\NCIC-2.lsm';
raw = tiffread(raw_pic_name);
green = raw.green;
height = raw.height;
width = raw.width;


seg = green > 0;
seg = imfill(seg,'holes');

b = uint16(edge(uint16(seg),'canny'));
%figure(1);
%imagesc(b);



L = bwlabel(b,8);
new = zeros(height,width);

%tag = L(15,129);%CIC-1;
%tag = L(72,11);%CIC-2;
%tag = L(18,102);%CIC-3;
%tag = L(128,104);%CIC-5;
%tag = L(389,357);%CIC-6;
%tag = L(237,413);%CIC-7;
%tag = L(350,217);%CIC-8;
%tag = L(6,182);%CIC-9;
%tag = L(362,227);%CIC-10;
%tag = L(306,148);%NCIC-1
tag = L(70,185);%NCIC-2
%tag = L(30,179);%NCIC-3
%tag = L(298,305);%NCIC-4
%tag = L(163,33);%NCIC-5

new1 = find(L==tag);
for k=1:1:size(new1,1)
    new(new1(k)) = 1;
end
%figure(2);
%imagesc(new);

boundary = find(new==1);
x = zeros(size(boundary,1),1);
y = zeros(size(boundary,1),1);
for k=1:1:size(boundary,1)
    if mod(boundary(k,1),height) ~= 0
        x(k,1) = fix(boundary(k,1)/height) + 1;
        y(k,1) = boundary(k,1) - (x(k,1)-1)*height;
    else
        y(k,1) = height;
        x(k,1) = boundary(k,1)/height;
    end
end

p = polyfit(x,y,1);
degree = atand(p(1));

disp(degree);

[x1,y1] = resort_boundary(x,y);
length = size(x1,1);

for k=200:1:299
    tempx(k-199,1) = x1(k,1);
    tempy(k-199,1) = y1(k,1);
end

new1 = zeros(height,width);
for k=1:1:100
    new1(tempy(k,1),tempx(k,1)) = 1;
end
%figure(3);
%imagesc(new1);



[degree1,mask] = chop_mask(tempx,tempy,height,width,100,1);

p.Nimages   = 1;  
p.yred      = 1:height;
p.xred      = 1:width; 
p = param_example(p); 


im3 = zeros(1,height,width);
im3(1,:,:) = green .* uint16(mask);
data = fire(p,im3,1);
%figure(4);
%imagesc(mask);
disp(degree1);

fiber = data.F;
points = data.X;
fiber_degree = zeros(size(fiber,2),1);
diff_degree = zeros(size(fiber,2),1);

for k=1:1:size(fiber,2)
    temp = fiber(k).v;
    tx = zeros(size(temp,2),1);
    ty = zeros(size(temp,2),1);
    for j=1:1:size(temp,2)
        tx(j) = points(temp(j),1);
        ty(j) = points(temp(j),2);
    end
    temp_p = polyfit(tx,ty,1);
    fiber_degree(k,1) = atand(temp_p(1));
    diff_degree(k,1) = abs(degree1 - fiber_degree(k,1));
    if diff_degree(k,1) > 90
        diff_degree(k,1) = 180-diff_degree(k,1);
    end
end
%}

%{
xlswrite('e:\Anju\Anju.xls',diff_degree,'Sheet2','P2');
xlswrite('e:\Anju\Anju.xls',fiber_degree,'Sheet3','P2');
m = mean(diff_degree);
s = std(diff_degree);
xlswrite('e:\Anju\Anju.xls',m,'Sheet4','P2');
xlswrite('e:\Anju\Anju.xls',s,'Sheet4','P3');
%}

figure(3);
hist(fiber_degree);



figure(4);
hist(diff_degree);

%[a,b,c,d] = EMSeg(fiber_degree,2);



%{
new = zeros(height,width);
for k=50:1:150
    new(y1(k,1),x1(k,1)) = 1;
end
figure(3);
imagesc(new);
%}
