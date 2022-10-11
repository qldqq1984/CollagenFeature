function [id,Num_T,Num_S] = getTSnum(fileName)

id = 0;
Num_T = 0;
Num_S = 0;
tag1 = strfind(fileName,'-');
if isempty(tag1) ~= 1
    id = str2double(fileName(1:tag1-1));
    tag_j = strfind(fileName,'j');
    tag_t = strfind(fileName,'t');
    tag_end = strfind(fileName,'.');
    if isempty(tag_j) == 1
        Num_T = Inf;
    elseif isempty(tag_t) == 1
        Num_S = Inf;
    else
        if tag_j < tag_t
            Num_S = tag_t-tag_j-1;
            Num_T = tag_end-tag_t-1;
        else
            Num_T = tag_j-tag_t-1;
            Num_S = tag_end-tag_j-1;
        end
    end
end


