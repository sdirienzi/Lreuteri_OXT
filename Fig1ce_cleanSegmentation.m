function cleanMask = cleanSegmentation(I, zscore, minArea, maxEcc)
    % CLEANSEGMENTATION cleans the segmentation using region properties. 
    % cleanMask = CLEANSEGMENTATION(I, zscore, minArea, maxEcc)
    % INPUT
    %   I: input mask, where 1 is cell (foreground). 
    %   zscore (optional): area z-score threshold, default: 2 (std). 
    %   minArea (optional): minimum area to keep, default: 30 (pixels).
    %   maxEcc (optional): maximum eccentricity to keep, default: 0.9.
    %       Circle has an eccentricity of 0. 
    % OUTPUT:
    %   cleanMask: cleaned logical mask. Foreground is true. 
    % 
    % Zhuohe Liu & Jihwan Lee
    % St-Pierre Lab, July 2019
    
    if nargin < 4
        maxEcc = 0.9; 
    end
    if nargin < 3
        minArea = 30;
    end
    if nargin < 2
        zscore = 2;
    end
    % convert files to black and white. 
    I = I == 1;      % 1 is cell == true

    I = imclearborder(I, 4);    % remove border ROIs
    
    % Eliminate too big or too small
    cc2 = bwconncomp(I, 4);     % we need to call this function to consider connected pixels. 
    stats0 = regionprops(cc2, 'Area', 'Eccentricity');
%     L = bwlabel(I, 4);        % gpuArray
%     stats0 = regionprops(L, 'Area', 'Eccentricity', 'PixelIdxList');  % gpuArray
    allArea = [stats0.Area];
    areaZscore = (allArea - mean(allArea))./std(allArea);
    unacceptedArea = abs(areaZscore) > zscore | ...
        allArea < minArea | ...
        [stats0.Eccentricity] > maxEcc;
%     sublist = cat(1, stats0(unacceptedArea).PixelIdxList);     % gpuArray
    sublist = cc2.PixelIdxList(unacceptedArea);
    sublist = cat(1, sublist{:}); 
    I(sublist) = false;
    cleanMask = I;     
end