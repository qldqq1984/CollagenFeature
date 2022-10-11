function [ori,stats] = getSHGfea_Orientation(mask,SHG,step)

%step = 200;
count = 0;
ori = [];
iso = [];
for i=1:step:size(mask,1)
    for j=1:step:size(mask,2)
        if i+step > size(mask,1) || j+step > size(mask,2)
            continue;
        end
        temp1 = mask(i:i+step-1,j:j+step-1);
        if (sum(sum(temp1))/(step*step)) < 0.8
            continue;
        else
            temp2 = SHG(i:i+step-1,j:j+step-1);
            if mean(mean(temp2)) < 20
                continue;
            end
            count = count + 1;
            %[ori(count),~,~,~,F,~]=ftmethod2(temp2);
            [ori(count),~,~,~,F,~]=ftmethod2(temp2);
            iso(count) = min(F(1,1),F(2,2))/max(F(1,1),F(2,2));
            %iso(count) = std(angle)/mean(angle);
            %if F(1,1)/F(2,2) > theR && F(1,1)/F(2,2) < 1/theR
            %    iso(count) = 1;
            %else
            %    iso(count) = 0;
            %end
        end
    end
end
stats = histogram_features(iso,10);

        


