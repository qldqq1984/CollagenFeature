function Cur = calc_fibercur(X,F,R)
%CALC_FIBERLEN - calculates the length of the fibers and returns them in
%the form of a field in F and a Length array

Cur = zeros(length(F),1);
for fi=1:length(F)
    fv = F(fi).v;
    temp_x = zeros(length(fv),1);
    temp_y = zeros(length(fv),1);
    for j=1:length(fv)
        v = fv(j);
        temp_x(j) = X(v,1);
        temp_y(j) = X(v,2);
    end
    Cur(fi,1) = get_curvature(temp_x,temp_y);
end
