function [Fea_T,Fea_S] = getSHGfea(channel1,tissue,Num_T,Num_S,MMP,ff)

Fea_T = [];
Fea_S = [];

[L,num] = bwlabel(tissue);
%{
SHG = channel1.*double(tissue);
ss = SHG(tissue == 1);
mss = mean(ss);
%}
fea = zeros(num,31);
Iseg = zeros(size(tissue));
for K=1:1:num
    mask = L == K;
    SHG = channel1.*double(mask);
    ss = SHG(mask == 1);
    mss = mean(ss);
    SHG = medfilt2(SHG,[3 3]);
    %collagen_mask = EMSeg(mat2gray(SHG),2);
    if mss < 4095*0.05
        collagen_mask = SHG > 4095*0.02;
    elseif mss < 4095*0.1
        collagen_mask = SHG > 4095*0.05;
    else
        collagen_mask = SHG > 4095*0.1;
    end
    %collagen_mask = double(otsu(SHG,2));
    %collagen_mask = imresize(collagen_mask,0.5);
    collagen_mask = bwareaopen(collagen_mask,10);
    Iseg = Iseg + collagen_mask;
    %imwrite(collagen_mask,fname1);
    
    [height, width] = size(collagen_mask);
    if sum(sum(collagen_mask))/sum(sum(mask)) < 0.001
        continue;
    else
        
        p.Nimages   = 1;  
        p.yred      = 1:height;
        p.xred      = 1:width; 
        p = param_example(p); 
        im3 = zeros(1,height,width);
        im3(1,:,:) = collagen_mask;
        seg = zeros(1,height,width);
        seg(1,:,:) = logical(collagen_mask);
        %data = fire(p,im3,2,seg);
        %M = network_stat(data.Xa,data.Fa,data.Va,data.Ra);
        %M.area = sum(sum(collagen_mask))/sum(sum(mask));
        %MMP = 146/512*2;
        %fea(K,1) = M.area;
        fea(K,1) = sum(sum(collagen_mask))/sum(sum(mask));
        %{
        fea(K,2) = M.fiber_num/(sum(sum(mask))*MMP*MMP);
        fea(K,3) = mean(M.L.*MMP);
        fea(K,4) = std(M.L.*MMP);
        fea(K,5) = mean(M.RF.*MMP);
        fea(K,6) = std(M.RF.*MMP);
        fea(K,7) = mean(M.fstr);
        fea(K,8) = std(M.fstr);
        fea(K,9) = M.xlinkdens;
        fea(K,10) = mean(M.xlinkspace);
        fea(K,11) = std(M.xlinkspace);
        %}
        %{
        fea(K,12) = sum(sum(mask));
        [~,iso] = getSHGfea_Orientation(mask,SHG,50);
        if isempty(iso) ~= 1
            %fea(K,13) = length(find(iso == 0));%/length(find(iso == 1));
            fea(K,13) = mean(iso);
            fea(K,14) = length(iso);
        end
        [~,iso] = getSHGfea_Orientation(mask,SHG,100);
        if isempty(iso) ~= 1
            %fea(K,14) = length(find(iso == 0));%/length(find(iso == 1));
            fea(K,15) = mean(iso);
            fea(K,16) = length(iso);
        end
        [~,iso] = getSHGfea_Orientation(mask,SHG,150);
        if isempty(iso) ~= 1
            %fea(K,16) = length(find(iso == 0));%/length(find(iso == 1));
            fea(K,17) = mean(iso);
            fea(K,18) = length(iso);
        end
        [~,iso] = getSHGfea_Orientation(mask,SHG,200);
        if isempty(iso) ~= 1
            %fea(K,17) = length(find(iso == 0));%/length(find(iso == 1));
            fea(K,19) = mean(iso);
            fea(K,20) = length(iso);
        end
        %}
        [~,stats] = getSHGfea_Orientation(mask,SHG,50);
        fea(K,2:7) = stats;
        [~,stats] = getSHGfea_Orientation(mask,SHG,100);
        fea(K,8:13) = stats;
        [~,stats] = getSHGfea_Orientation(mask,SHG,150);
        fea(K,14:19) = stats;
        [~,stats] = getSHGfea_Orientation(mask,SHG,200);
        fea(K,20:25) = stats;
        [~,stats] = getSHGfea_Orientation(mask,SHG,250);
        fea(K,26:31) = stats;
    end
    
end

if Num_T == Inf 
    Fea_T = fea;
elseif Num_S == Inf
    Fea_S = fea;
else
    [~,list] = sort(fea(:,1),'descend');
    Fea_S = fea(list(1:Num_S),:);
    Fea_T = fea(list(Num_S+1:end),:);
end
%imwrite(Iseg,ff);    
    




