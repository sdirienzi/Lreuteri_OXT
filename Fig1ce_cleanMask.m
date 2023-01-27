    function I = cleanMask(mask)
    %mask = imread(path);
    mask = mask == 1;
    %crop_region = zeros(2044,549);
    %new_mask = cat(2, mask(:,1:1499), crop_region); 
    new_mask = mask; 
    %clean_mask = cleanSegmentation(mask, 100, 500, 0.95); 

%     %Apply imimposmin function 
%     CC =bwconncomp(new_mask);
%     stats = regionprops(CC, 'Area');
%     area = [stats.Area];
%     fusedCandidates = area > mean(area) + 2.0 * std(area);
%     sublist = CC.PixelIdxList(fusedCandidates);
%     sublist = cat(1,sublist{:}); 
%     fusedMask = false(size(new_mask));
%     fusedMask(sublist) = 1; 
%     %imshow(fusedMask); 
% 
%     s = round(1.2*sqrt(mean(area))/pi);
%     nucmin = imerode(fusedMask, strel('disk',s)); 
%     %imshow(nucmin); 
% 
%     outside= ~imdilate (fusedMask, strel('disk',1)); 
%     %imshow(outside); 
% 
%     basin = imcomplement(bwdist(outside));
%     basin = imimposemin(basin, nucmin | outside);
%     %pcolor(basin); shading flat; 
%     L = watershed(basin);
%     rgb = label2rgb(L,'jet',[.5 .5 .5]);
%     %imshow(rgb); 
%     combinedMask = L > 1 | (new_mask - fusedMask); 
%     %imshow(combinedMask); 

    %I = cleanSegmentation(new_mask, 10000, 500, 1.0); 
    %I = bwlabel(I);
    I = new_mask;
    
    mask = bwlabel(mask);
    %figure(); imshow(label2rgb(mask,'jet', 'k', 'shuffle'));
    %figure(); imshow(label2rgb(I,'jet', 'k', 'shuffle'));
    
    %MaskArray{count} = clean_mask;
    %count = count + 1 
    
