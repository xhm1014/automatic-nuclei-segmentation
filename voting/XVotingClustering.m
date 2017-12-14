function [s2]= XVotingClustering(bw_i,img,Para)
% Input
%     -bw_i the binary image
%     -img the original image for observation
% Output
%     -RC_i the mask image with markers
% Program written by Hongming Xu
% Electrical and computer engineering, university of alberta

%% parameters setting

rmin=Para.rmin;
rmax=Para.rmax;
delta=Para.delta;
bandwidth=Para.bandwidth;
Sigma=Para.Sigma;


%% compute the number of connected components
[m,n] = size(bw_i);
c1 = bwconncomp(bw_i);
for i = 1:c1.NumObjects
    bw_cur = bw_i;
    bw_cur(:) = 0;
    bw_cur(c1.PixelIdxList{i}) = 1;
    Mul_bw{i} = bw_cur;
end
num1=length(Mul_bw);

%% find the high gradients points
ltg = cell(1,num1);    % the high gradients
lindex = cell(1,num1); % the positions of high gradients
dxy=cell(1,num1);

for t = 1:length(Mul_bw)
    bw_cur = im2double(Mul_bw{t});
    [sIx,sIy] = gradient(bw_cur);
    g = sqrt(sIx.^2+sIy.^2);
    ind = find(g);
    
    %% fast version
    ltg{t}=g(ind);
    [r,c]=ind2sub([m,n],ind);
    lindex{t}=[r';c'];
    
    dx = sIx(ind)./g(ind);
    dy = sIy(ind)./g(ind);
    dxy{t}=[dx';dy'];
    %     tg = min(g(ind));  % the threshold to find boundary
    %     for i = 1:2:m
    %         for j = 1:2:n
    %             if g(i,j) >= tg
    %                ltg{t} = [ltg{t},g(i,j)];
    %                lindex{t} = [lindex{t},[i;j]];
    %             end
    %         end
    %     end
end

%% ellipse algorithm
PAR = zeros(length(lindex),5);
for i = 1:length(lindex)
    X = lindex{i}(2,:);
    Y = lindex{i}(1,:);
    PAR(i,:) = fitellipse(X,Y);
end
[x1,y1] = meshgrid(1:n,1:m);
Cen = [;];
num = [];
for i = 1:size(PAR,1)
%    t = linspace(0,pi*2);
%    x = PAR(i,3) * cos(t);
%    y = PAR(i,4) * sin(t);
%    nx = x*cos(PAR(i,5))-y*sin(PAR(i,5)) + PAR(i,1);
%    ny = x*sin(PAR(i,5))+y*cos(PAR(i,5)) + PAR(i,2);
    
    
    K_bw = zeros(m,n);
    X = (x1(:) - repmat(PAR(i,1),m*n,1)).*cos(PAR(i,5)) +(y1(:) - repmat(PAR(i,2),m*n,1)).*sin(PAR(i,5));
    Y = -(x1(:) - repmat(PAR(i,1),m*n,1)).*sin(PAR(i,5)) + (y1(:) - repmat(PAR(i,2),m*n,1)).*cos(PAR(i,5));
    ind = find(((X.^2)./(PAR(i,3).^2)+(Y.^2)./(PAR(i,4).^2)) <= 1);
    K_bw(ind) = 1;
    
    t_bw = K_bw & Mul_bw{i};
    e1 = double(sum(t_bw(:)))/double(sum(K_bw(:)));
    
    r_bw = bitxor(K_bw,Mul_bw{i});
    e2 = double(sum(r_bw(:)))/double(sum(K_bw(:)));
    %    s1 = regionprops(Mul_bw{i},'Area');
    
    
    thr1=0.91;thr2=0.18;
    
    if e1 > thr1 && e2 < thr2
        num = [num,i];
        s = regionprops(Mul_bw{i},'centroid');
        Cen = [Cen;s.Centroid];
    end
end
for i = length(num):-1:1
    Mul_bw(:,num(i)) = [];
    lindex(:,num(i)) = [];
    ltg(:,num(i)) = [];
    dxy(:,num(i)) = [];
end


%% coarse erosion
bw_e = XCondition_Erosion(Mul_bw);
%bw_e=Mul_bw;

%% the new voting algorithm
J = zeros(m,n); % the voting image
for t = 1:length(Mul_bw)
    for i = 1:length(ltg{t})
        r = lindex{t}(1,i);
        v = lindex{t}(2,i);
        
        dx=dxy{t}(1,i);dy=dxy{t}(2,i);
        
        % compute shift Gaussian kernel
        ux = v+(rmax+rmin)*(dx)/2;
        uy = r+(rmax+rmin)*(dy)/2;
        
        [rt,ct] = find(bw_e{t});
        TX = [rt ct];
        Ty = [r v];
        dis = sqrt(sum(abs(TX - repmat(Ty,length(rt),1)).^2,2));
        tindex = find(dis>rmin & dis<rmax);
        
        TX2 = TX(tindex,:); % points in the fan
        TV = TX2 - repmat(Ty,size(TX2,1),1); % points within the angle
        t_mat = TV(:,1);
        TV(:,1) = TV(:,2);
        TV(:,2) = t_mat;
        Td = repmat([dx dy],size(TV,1),1);
        tc = sum(Td.*TV,2)./(sqrt(sum(Td.*Td,2)).*sqrt(sum(TV.*TV,2)));
        c_delta = acos(tc);
        cindex = find(c_delta < delta/2);
        findex = TX2(cindex,:);
        
        T_u = repmat([uy ux],size(findex,1),1);
        gau = 1/(sqrt(2*pi)*Sigma).*exp(-(sum((findex-T_u).^2,2))./(2*Sigma^2));
        ind = sub2ind([m n],findex(:,1)',findex(:,2)');
        J(ind) = J(ind)+15*ltg{t}(i)*gau';
        %         disp(i);
    end
end
%   show(J);

%% mean shift clustering
% the coordinates are the clustering feature
X = [;];
[row,col] = find(J >= 0.5);
for k = 1:length(row)
    i = row(k);
    j = col(k);
    temp = round(J(i,j));
    xtemp = repmat([j;i],1,temp);
    X = [X,xtemp];
end


% mean shift algorithm
[clustCent,point2cluster,clustMembsCell] = MeanShiftCluster(X,bandwidth);
s2=[];
if length(Cen)>0;
    s2=[s2;round(Cen(:,2)),round(Cen(:,1))];
end
if length(clustCent)>0
    s2=[s2;round(clustCent(2,:))',round(clustCent(1,:))'];
end
%% show results
% figure(1),imshow(img);
% hold on,
% if length(Cen)>0
%       plot(Cen(:,1),Cen(:,2),'gx','MarkerSize',7,'LineWidth',2);         % plot the center by elliptical model
%
% end
% if length(clustCent)>0
%      plot(clustCent(1,:),clustCent(2,:),'r+','MarkerSize',7,'LineWidth',2); % plot the center by mean shift
% end
%  hold off;
%
%
%
%
% RC_i = zeros(m,n);
% if length(clustCent)>0
%     linearInd = sub2ind([m n],round(clustCent(2,:)),round(clustCent(1,:)));
%     RC_i(linearInd) = 1;
% end
% if length(Cen)>0
%     linearInd = sub2ind([m n],round(Cen(:,2)),round(Cen(:,1)));
%     RC_i(linearInd) = 1;
% end


