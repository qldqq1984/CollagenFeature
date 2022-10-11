function [degree,mask] = chop_mask(x,y,height,width,depth,method)

mask = zeros(height,width);


p = polyfit(x,y,1);
degree = atand(p(1));

if degree == 0
    degree1 = 90;
elseif degree > 0 
    degree1 = degree - 90;
else
    degree1 = degree + 90;
end

length = size(x,1);


kk = tand(degree1);
bb1 = y(length,1) - kk*x(length,1);
bb2 = y(1,1) - kk*x(1,1);

distance = round(depth*sind(abs(degree)));
%disp(distance);
%pause();

if method == 1
    tt = 1;
    for i=x(1,1):-1:x(1,1)-distance
        xx(tt,1) = i;
        yy(tt,1) = round(i*kk + bb2);
        tt = tt+1;
    end
    
    sx = xx(distance,1);
    sy = yy(distance,1);
    clear xx yy;
    
    tt = 1;
    for i=x(length,1):-1:x(length,1)-distance
        xx(tt,1) = i;
        yy(tt,1) = round(i*kk + bb1);
        tt = tt+1;
    end
    
    ex = xx(distance,1);
    ey = yy(distance,1);
    clear xx yy;

    c = [x(1,1) x(length,1) ex sx];
    r = [y(1,1) y(length,1) ey sy];
    mask = roipoly(mask,c,r);
    
elseif method == 2
    tt = 1;
    for i=x(1,1):1:x(1,1)+distance
        xx(tt,1) = i;
        yy(tt,1) = round(i*kk + bb2);
        tt = tt+1;
    end
    
    sx = xx(distance,1);
    sy = yy(distance,1);
    clear xx yy;
    tt = 1;
    for i=x(length,1):1:x(length,1)+distance
        xx(tt,1) = i;
        yy(tt,1) = round(i*kk + bb1);
        tt = tt+1;
    end
    
    ex = xx(distance,1);
    ey = yy(distance,1);
    clear xx yy;

    c = [x(1,1) x(length,1) ex sx];
    r = [y(1,1) y(length,1) ey xy];
    mask = roipoly(mask,c,r);
end

