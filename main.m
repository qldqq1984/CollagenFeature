addpath(genpath(pwd));

clear all
close all 

DataSets_dir=[pwd,'\photo\*.tif',];
S = dir(DataSets_dir);
filearray = dir(DataSets_dir);
fea = zeros(length(filearray),8);
Fname = cell(length(filearray),1);

for K=1:1:length(filearray)

    fileName = [pwd,'\photo\' filearray(K).name];
    fileName1 = [pwd,'\Mask\' filearray(K).name(1:end-3) 'tif'];
    disp(filearray(K).name(1:end-4));
    Fname{K,1} = filearray(K).name(1:end-4);
  
    data=imread(fileName);
     SHG=data(:,:,2);
     channel1 =SHG;
     
    threshold = 45;   %%threshold (0-255)

    if max(max(channel1)) > 255
        channel1 = uint8(channel1./2^8);
        channel2 = medfilt2(channel1,[3 3]);    
        collagen_mask = channel2 > threshold;  
        %end
    else
        channel1=uint8(channel1./2^0);
        channel2 = medfilt2(channel1,[3 3]);
        collagen_mask = channel2 > threshold;

    end
    collagen_mask = bwareaopen(collagen_mask,5);
    imwrite(collagen_mask,fileName1);
    height = size(collagen_mask,1);
    width = size(collagen_mask,2);

    if (sum(sum(collagen_mask))/(height*width)) < 0.01
        fea(K,1) = sum(sum(collagen_mask))/(height*width);
    else
%     %morphology features
        p.Nimages   = 1;  
        p.yred      = 1:height;
        p.xred      = 1:width; 
        p = param_example(p); 
        im3 = zeros(1,height,width);
        im3(1,:,:) = collagen_mask;
        seg = zeros(1,height,width);
        seg(1,:,:) = logical(collagen_mask);
        data = fire(p,im3,2,seg);
        M = network_stat(data.Xa,data.Fa,data.Va,data.Ra);

        MMP =  0.263; %%%Resolution ratio

        fea(K,1) = sum(sum(collagen_mask))/(height*width);
        fea(K,2) = M.fiber_num/((height*width)*MMP*MMP);
        fea(K,3) = mean(M.L.*MMP);
        fea(K,4) = mean(M.RF.*MMP);
        fea(K,5) = mean(M.fstr);
        fea(K,6) = M.xlinkdens;
        fea(K,7) = mean(M.xlinkspace);
        [ori,~,~,~,F,~]=ftmethod2(channel1);
        fea(K,8) = min(F(1,1),F(2,2))/max(F(1,1),F(2,2));

      
    end
  
end


