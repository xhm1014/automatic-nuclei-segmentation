function [bw_e] = XCondition_Erosion(bw)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Input
%     -bw the binary images
% Output
%    -bw_e the binary image after condition erosion
% Program written by Hongming Xu

%% build the four coarse erosion structures 
% strel_c1 = strel([0 0 0 0 0 0 0;0 0 1 1 1 0 0;...
%     0 1 1 1 1 1 0; 1 1 1 1 1 1 1; 0 1 1 1 1 1 0;0 0 1 1 1 0 0;0 0 0 0 0 0 0]);
% strel_c2 = strel([0 0 0 0 0 0 0;0 0 1 1 1 0 0;...
%     0 1 1 1 1 1 0; 1 1 1 1 1 1 1; 0 1 1 1 1 1 0;0 0 1 1 1 0 0;0 0 0 0 0 0 0]');
% 
% strel_c3 = strel([0 0 1 1 1 0 0;...
%     0 1 1 1 1 1 0; 1 1 1 1 1 1 1; 1 1 1 1 1 1 1;1 1 1 1 1 1 1;0 1 1 1 1 1 0;0 0 1 1 1 0 0]);
strel_c6 = strel('disk',3,0);

%% perform the coarse erosion    
for i = 1:length(bw)
    bw_cur = bw{i};
    while sum(bw_cur(:)) > round(sum(bw{i}(:))*0.5)
        c1 = bwconncomp(bw_cur);
        if c1.NumObjects == 1
%                 bw_cur = imerode(bw_cur,strel_c1);
%                 bw_cur = imerode(bw_cur,strel_c2);
                bw_cur = imerode(bw_cur,strel_c6);
        else
            bw_cur(:) = 0;
            for j = 1:c1.NumObjects
                bw_j = bw_cur;
                bw_j(:) = 0;
              
                bw_j(c1.PixelIdxList{j}) = 1;
                if sum(bw_j(:)) > 100
%                         bw_j = imerode(bw_j,strel_c1);
%                         bw_j = imerode(bw_j,strel_c2);
                        bw_j = imerode(bw_j,strel_c6);
                end
                bw_cur = bw_cur + bw_j;
            end
        end
    end
    bw{i} = bw_cur;
end

bw_e = bw;
