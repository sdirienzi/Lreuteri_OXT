%% ReadMe
% Inputs: 
% (1) Raw nucleus stained DAPI image
% (2) Raw plasma membrane stained image
% (3) Coressponding binary nucleus mask -> Generated with Ilastik; label cells with value 1 and background with value 2. 
% (4) Pixel Threhold value to generate plasma membrane mask. -> to detect
% epithelial cells
% (5) Visual cue -> if visuals are necessary make value = 1 if not = 0. 
% Outputs:
% (1) Total number of nucleus in the image 
% Dependency:
% (1) cleanSegmentation.m
% (2) cleanMask.m 
% -> Have these two codes in the working directory. 
% This code is written by Jihwan Lee, May 2021.
% Warning : 
% (1) None
%% 
function [nucleus_count_by_area] = countCell_usingArea (raw_dapi_image, raw_membrane_image, membrane_mask_threshold, nucleus_mask_threshold, average_nuc_size, visual_cue)
    dapi_image = imread(raw_dapi_image);
    cy5_image = imread(raw_membrane_image);
    %dapi_image_mask = imread(mask_image);
    
% Generate masks 
% Define region that corresponds to epithelial cells
    cy5_mask = cy5_image > membrane_mask_threshold;
    se = strel('disk', 5, 0);
    cy5_mask_dilate_fill = imdilate(cy5_mask, se);
    %imshow(cy5_mask_dilate , []);
    %cy5_mask_dilate_fill = imfill(cy5_mask_dilate, 'holes');

% Generate binary mask for nucleus
    %dapi_image_mask_clean =cleanMask(dapi_image_mask);
    dapi_image_mask_clean = dapi_image > nucleus_mask_threshold;
% Locate where nuc mask and epi mask overalp 
nuc_epi_comb_mask = cy5_mask_dilate_fill & dapi_image_mask_clean;
nuc_epi_comb_mask_clean = cleanMask(nuc_epi_comb_mask);
CC =bwconncomp(nuc_epi_comb_mask_clean);
stats_nuc_epi = regionprops(CC,'Centroid', 'Area');
total_nuc_area = sum([stats_nuc_epi.Area]);
%average_nuc_size = 3225 % in pixel
nucleus_count_by_area = round(total_nuc_area ./ average_nuc_size); 
    
% %{ Count total cells   
%     
%     CC =bwconncomp(dapi_image_mask_clean);
%     stats1 = regionprops(CC,'Centroid', 'Area'); 
%     centroids = cat(1,stats1.Centroid);
%     
%     total_count = size(stats1, 1); 
    
% % Count epithelial cells  
%     round_centroids = round(centroids);
% 
%     
%     index_list = [];
%     for i = 1:size(round_centroids,1)
%         if cy5_mask_dilate_fill(round_centroids(i,2), round_centroids(i,1)) == 1;
%         index_list =[index_list; i];
%         else 
%         continue 
%         
%         end 
%     end 
%     new_stats = stats1(index_list);
%     
%    
%     epi_count = size(new_stats, 1); 
    

% Count big cells

%figure(); histogram([new_stats.Area] ./ mean([new_stats.Area]))
%big_cell = find([new_stats.Area] ./ mean([new_stats.Area]) > 2.0);
%big_cell_stats = new_stats(big_cell);
%big_epi_count = size(big_cell_stats, 1); 
%big_epi_count_corrected = sum(round([big_cell_stats.Area] ./ mean([new_stats.Area])));
%if big_epi_count ./ epi_count > 0.1 
    %warning('Too many big nucleus. Check segmentation mask')
%end 
%epi_count_corrected = epi_count - big_epi_count + big_epi_count_corrected;
    

% Visualize
if visual_cue == 1 
    
    figure(); 
    subplot (2,3,1)
    imshow(dapi_image , []);
    title('Nucleus staining');

    subplot(2,3,2)
    imshow(cy5_image , []);
    title('Plasma membrane image');

    subplot(2,3,3)
    imshow(cy5_mask_dilate_fill , []);
    title('Epithelial cell mask');
    

    subplot(2,3,4)
    
    imshow(dapi_image_mask_clean, []);
    title('Nuclei mask');
    hold on
    
   
    subplot(2,3,5)
    imshow(nuc_epi_comb_mask_clean, []);
    
    title('Nuclei + epi combined mask')
elseif visual_cue == 0
    
end 







    