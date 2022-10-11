clear all;
close all;
clc;

addpath(genpath('C:\Shuoyu\code\'));
%{
raw = tiffread('d:\Shuoyu\TS3.lsm');
height = raw.height;
width = raw.width;

mask = imread('d:\Shuoyu\TS3.tif');
mask = mask./255;
%}
I = imread('F:\Finish\Finish new\Week6\311007-16-1.mdb\17_EM_combine2.tif');
mask = I(2001:2500,2001:2500);
mask = mask./255;
%imagesc(a)
[height width] = size(mask);

p.Nimages   = 1;  
p.yred      = 1:height;
p.xred      = 1:width; 
p = param_example(p); 
im3 = zeros(1,height,width);
im3(1,:,:) = mask;
seg = zeros(1,height,width);
seg(1,:,:) = logical(mask);
data = fire(p,im3,2,seg);
M = network_stat(data.Xa,data.Fa,data.Va,data.Ra);
figure(3);
imagesc(mask);

length = size(data.Fa,2);
count = 0;
count1 = 0;
for k=1:1:length
    if isempty(data.Fa(k).f) == 1
        %disp(k);
        count = count + 1;
        dist(count) = 1;
    else
        count1 = count1 + 1;
        agg(count1) = 1;
    end
end
a = find(dist == 1);
data_d.Xa = data.Xa;
data_d.Fa = data.Fa(a');
data_d.Va = data.Va;
data_d.Ra = data.Ra;

b = find(agg == 1);
data_a.Xa = data.Xa;
data_a.Fa = data.Fa(b');
data_a.Va = data.Va;
data_a.Ra = data.Ra;

M_d = network_stat(data_d.Xa,data_d.Fa,data_d.Va,data_d.Ra);
M_a = network_stat(data_a.Xa,data_a.Fa,data_a.Va,data_a.Ra);
%disp(count);



