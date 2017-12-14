%%----------------------%%
% Input
%     R: image channel (R,G,B,or Gray channel)
%     fm: the foregrond mask (generally exclued while pixels)
%     ac: remove non-nuclei regions (i.e., noisy regions)
% Output:
%     Nmask: the foreground is nuclei which clump together or have
%    concave shape
%    cs4,rs4: the centers of isolated nuclei whith convex shape
%    A3: the isolated nuclei with convex shape
% Written by Hongming Xu
% ECE U of A
% feel free to use it
%%---------------------%%


function [Nmask]=XNucleiBinarization(R,fm,ac)

%% HGMR module
R_hat=255-R;
% Opening by reconstruction
% S = [0 0 1 1 1 0 0;...
%     0 1 1 1 1 1 0;...
%     1 1 1 1 1 1 1;...
%     1 1 1 1 1 1 1;...
%     1 1 1 1 1 1 1;...
%     0 1 1 1 1 1 0;...
%     0 0 1 1 1 0 0];

S = [0 0 1 0 0;...         %% this parameter should be adaptive changed according to nuclei size shape
     0 1 1 1 0;...
     1 1 1 1 1 ;...
     0 1 1 1 0;...
     0 0 1 0 0];
Re=imerode(R_hat,S);
fRe=imreconstruct(Re,R_hat);

% Closing-by-Reconstruction
fRerc=imcomplement(fRe);
fRerce=imerode(fRerc,S);
fRercbr=imcomplement(imreconstruct(fRerce,fRerc));

R=255-fRercbr;
%R=LmakeHomoonRC(R);
%% Global thresholding
temp=R(fm);
%TR=graythresh(temp);
%TR=0.6;                 %% this value is more robust than Otus's method Dataset II
TR=0.67;   %%DatasetI
Rlogical=im2bw(R,TR);
RClogical=Rlogical|(~fm); % obtain a binary image with nuclei as the foreground
C_mask1=imopen(~RClogical,strel('disk',2));
C_mask1=bwareaopen(C_mask1,ac,8);
%C_mask1=imfill(C_mask1,'holes');

%% local threshold segmentation
[label2,n2]=bwlabel(C_mask1);
stats2=regionprops(label2,'BoundingBox');
C_mask2=C_mask1;
%imshow(R)
for j=1:n2
    x=floor(stats2(j).BoundingBox(1));
    y=floor(stats2(j).BoundingBox(2));
    w=floor(stats2(j).BoundingBox(3));
    h=floor(stats2(j).BoundingBox(4));
    if y<1
        y=1;
    end
    if y+h>size(R,1)
        y2=size(R,1);
    else
        y2=y+h;
    end
    if x<1
        x=1;
    end
    if x+w>size(R,2)
        x2=size(R,2);
    else
        x2=x+w;
    end
    tr=R(y:y2,x:x2);
    rr=im2bw(tr,graythresh(tr));   %% Otus's method
    C_mask2(y:y2,x:x2)=C_mask1(y:y2,x:x2)&(~rr);
%      hold on,                    %% Drawing bounding boxes
%      plot(stats2(j).BoundingBox(1),stats2(j).BoundingBox(2),'r*');
%      rectangle('Position',[stats2(j).BoundingBox(1),stats2(j).BoundingBox(2),stats2(j).BoundingBox(3),stats2(j).BoundingBox(4)],...
%         'EdgeColor','g','LineWidth',2);
end
C_mask2=imopen(C_mask2,strel('disk',2));
C_mask2=bwareaopen(C_mask2,round(ac*(2/3)),8);
Nmask=C_mask2;
%Nmask=imfill(C_mask2,'holes');

end