function [ROIs] = getROIs(locList,ROIrad,imSize)
%GETROIS get ROIs from list of localizations
%   Detailed explanation goes here

% build ROIs
ROIcent = round(locList);
%in main from LAIA, ROIcent is given as [X,Y]
nLoc   = size(ROIcent,1);
% imSize = [size(data,1),size(data,2)];

% storing of ROIs
% ROIS
%ROI is given as [x1,x2,y1,y2]
ROIs = zeros(nLoc,4);
ROIs(:,1) = round(ROIcent(:,1))-ROIrad;
ROIs(:,2) = round(ROIcent(:,1))+ROIrad;
ROIs(:,3) = round(ROIcent(:,2))-ROIrad;
ROIs(:,4) = round(ROIcent(:,2))+ROIrad;

%fix in case
if ROIs(:,1) <= 0
    ROIs(:,2) = ROIs(:,2)-ROIs(:,1)+1;
    ROIs(:,1) =1;
end
if ROIs(:,2) > imSize(2)
    
    ROIs(:,1) = ROIs(:,1)-(ROIs(:,2)-imSize(2));
    ROIs(:,2) =imSize(2);
end



if ROIs(:,3) <= 0
    ROIs(:,4) = ROIs(:,4)-ROIs(:,3)+1;
    ROIs(:,3) =1;
end

if ROIs(:,4) > imSize(1)
    
    ROIs(:,3) = ROIs(:,3)-(ROIs(:,4)-imSize(1));
    ROIs(:,4) =imSize(1);
end

%imLims receive imSize(2) first so it is also comparing XY with XY)
imLims = repmat([1 imSize(2) 1 imSize(1)],nLoc,1);
test = [ROIs(:,1)>=imLims(:,1), ROIs(:,2)<=imLims(:,2),...
        ROIs(:,3)>=imLims(:,3), ROIs(:,4)<=imLims(:,4)];
    
assert(all(test(:)), 'Problems, unexpected issues during ROI definition')
assert((ROIs(:,2)-ROIs(:,1))== ((ROIs(:,4)-ROIs(:,3))),'ROI not symmetric')
end

