function[M] = network_stat(X,F,V,R)
%NETWORK_STAT - measures properties of the network

%construct an edge matrix
    E = zeros(length(F),2);
    for fi=1:length(F)
        E(fi,:) = [F(fi).v(1) F(fi).v(end)];
    end

%count number of fibers
    N = length(E);
    M.fiber_num = N;
  
%calculate length of all the fibers    
    [F, L, RF] = calc_fiberlen(X,F,R);  
    M.L       = L;      
    M.avgL    = mean(L);   
    M.totL    = sum(L);
    M.RF = RF;
    M.avgR = mean(RF);
    M.totR = sum(RF);
    %CU = calc_fibercur(X,F,R);
    %M.CU = CU;
    %M.avgC = mean(CU);

%calculate  strightness = (length of the straight line connecting the fiber start point and the end point )/ (fiber length)    
    fnum = N;  % nuber of the extracted fibers
    fsp = zeros(fnum,1); % fiber start point
    fep = zeros(fnum,1); % fiber end point
    dse = zeros(fnum,1); % distance between start point and end point
    for i = 1:fnum
        fsp(i) = F(i).v(1);
        fep(i) = F(i).v(end);
        dse(i) = norm(X(fep(i),:)-X(fsp(i),:));
    end
    M.fstr = dse./M.L;
   
%{
%calculate volume fraction and density of network
    %M.Ldens   = M.totL/prod(max(X)-min(X));
    M.Ldens = 0;
    M.fibervol= sum(L.*pi.*RF.^2);
    M.vol     = prod(max(X)-min(X));
    %M.volfrac = M.fibervol/M.vol;
    M.volfrac = 0;
    M.dens    = M.volfrac*1360;
    1;
%}    
%{
%calculate angle orientation of fibers
    v1 = E(:,1);
    v2 = E(:,2);
    x1 = X(v1,:);
    x2 = X(v2,:); 
    
    M.angle_xz   = atan( (x2(:,3)-x1(:,3))./(x2(:,1)-x1(:,1)+eps) );
    M.angle_xy   = atan( (x2(:,2)-x1(:,2))./(x2(:,1)-x1(:,1)+eps) );
%}    

%calculate coordination number
    E = fiber2edge(F,V);
    coord = zeros(max(E(:)),1);
    for k=1:size(E,1)
        v = E(k,:);
        coord(v) = coord(v)+1;
    end
    ind = find(coord~=0);
    M.avgcoord = mean(coord~=0);
    M.coord    = coord;

%cross-link density and number
    xlink = 0;
    vflag = zeros(length(V),1);
    for i=1:length(V)
        if length(V(i).f)>1
            xlink = xlink+1;
            vflag(i) = 1;
        end
    end
    M.xlinkdens = xlink/M.totL;
    
%mean cross-link spacing distribution
    sp = 0;
    xlinkspace = zeros(xlink,1);
    
    for fi=1:length(F)
        v      = F(fi).v;
        ind    = find(vflag(v)==1);
        if length(ind) > 1 %there is more than one cross link in fiber
            for j=1:length(ind)-1
                sp = sp+1;
                i1 = ind(j);
                i2 = ind(j+1);
                len= 0;
                for ii=i1:i2-1
                    v1 = v(ii);
                    v2 = v(ii+1);
                    len = len+norm(X(v2,:)-X(v1,:));
                end
                xlinkspace(sp) = len;
            end
        end
    end
    if sp>0
        M.xlinkspace = xlinkspace(1:sp);
    else
        M.xlinkspace = []; %no cross-links in network
    end
    M.avexlinkspace = mean(M.xlinkspace);
%persistence length
%    [Lp3] = calc_persistent2(X,F);
%    M.Lp3 = Lp3;
%    1;
