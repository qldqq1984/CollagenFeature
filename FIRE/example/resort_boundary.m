function [x1,y1] = resort_boundary(x,y)

length = size(x,1);

temp_x = zeros(length,1);
temp_y = zeros(length,1);
tag = zeros(length,1);

temp_x(1,1) = x(1,1);
temp_y(1,1) = y(1,1);
tag(1,1) = 1;

for k=2:1:length
    find_1 = 0;
    for j=1:1:length
        if tag(j,1)~=1
            if (abs(x(j,1)-temp_x(k-1,1)) + abs(y(j,1)-temp_y(k-1,1))) == 1
                temp_x(k,1) = x(j,1);
                temp_y(k,1) = y(j,1);
                tag(j,1) = 1;
                find_1 = 1;
                break;
            end
        end
    end
    if find_1 == 1
        continue;
    else
        %find_2 = 0;
        for i=1:1:length
            if tag(i,1)~=1
                if (abs(x(i,1)-temp_x(k-1,1)) + abs(y(i,1)-temp_y(k-1,1))) == 2
                    temp_x(k,1) = x(i,1);
                    temp_y(k,1) = y(i,1);
                    tag(i,1) = 1;
                    %find_2 = 1;
                    break;
                end
            end
        end
    end
end
x1 = temp_x;
y1 = temp_y;
