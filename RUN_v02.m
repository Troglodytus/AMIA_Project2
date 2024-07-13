% RUN.m - Main script to classify different cell types on microscopy images

% Clear workspace and close all figures
clear; close all; clc;

% Set random seed for reproducibility
rng(42);

% Define the main directory containing images of different cell types
mainDir = 'C:\Users\Daniel\Documents\FHKärnten\Medical Engineering\2. Semester\Applied Medical Image Analysis\Project_Part2\bloodcells_4classes_500_test';
cellTypes = {'basophil', 'eosinophil', 'erythroblast', 'ig'};
numCellTypes = length(cellTypes);

% Define image size to resize to
resizeSize = [256, 256];

% Loop through each cell type folder
for i = 1:numCellTypes
    cellType = cellTypes{i};
    imageDir = fullfile(mainDir, cellType);
    imageFiles = dir(fullfile(imageDir, '*.jpg'));
    
    % Load one representative image from the folder
    imgPath = fullfile(imageDir, imageFiles(1).name);
    img = imread(imgPath);
    
    % Resize image to the predefined size
    imgResized = imresize(img, resizeSize);
    imgHSV = rgb2hsv(imgResized);

    % Convert image to grayscale for visualization
    imgGray = rgb2gray(imgResized);

    % Define the seed point (center of the image)
    seedPoint = [round(resizeSize(1)/2), round(resizeSize(2)/2)];
    
    % Region growing parameters
    hueLowThreshold = 0.3;
    hueHighThreshold = 0.99;
    brightnessLowThreshold = 0.25;
    brightnessHighThreshold = 0.85;
    saturationThreshold = 0.1;
    hueBackgroundThreshold = 0.5;
    maxHsvDistance = 0.1; % Maximum allowable HSV distance for new points
    
    % Initialize region growing
    regionMask = false(resizeSize);
    regionMask(seedPoint(1), seedPoint(2)) = true;
    growingList = seedPoint;
    hsvValues = squeeze(imgHSV(seedPoint(1), seedPoint(2), :))';
    checklist = hsvValues; % Initialize checklist with the seed point HSV values

    % Perform region growing for the initial 10 steps
    for step = 1:10
        if isempty(growingList)
            break;
        end
        
        currentPoint = growingList(1, :);
        growingList(1, :) = [];

        for k = -1:1
            for l = -1:1
                newRow = currentPoint(1) + k;
                newCol = currentPoint(2) + l;

                if newRow > 0 && newRow <= resizeSize(1) && newCol > 0 && newCol <= resizeSize(2)
                    if ~regionMask(newRow, newCol)
                        hueValue = imgHSV(newRow, newCol, 1);
                        saturationValue = imgHSV(newRow, newCol, 2);
                        brightnessValue = imgHSV(newRow, newCol, 3);
                        
                        % Add point to checklist regardless of criteria
                        checklist = [checklist; hueValue, saturationValue, brightnessValue];
                        
                        % Check if the point fits the cell criteria
                        if hueValue >= hueLowThreshold && hueValue <= hueHighThreshold && ...
                           brightnessValue >= brightnessLowThreshold && brightnessValue <= brightnessHighThreshold && ...
                           saturationValue > saturationThreshold
                            regionMask(newRow, newCol) = true;
                            growingList = [growingList; newRow, newCol];
                            hsvValues = [hsvValues; hueValue, saturationValue, brightnessValue];
                        end
                    end
                end
            end
        end
    end

    % Perform the main region growing with additional similarity check
    while ~isempty(growingList)
        currentPoint = growingList(1, :);
        growingList(1, :) = [];

        for k = -1:1
            for l = -1:1
                newRow = currentPoint(1) + k;
                newCol = currentPoint(2) + l;

                if newRow > 0 && newRow <= resizeSize(1) && newCol > 0 && newCol <= resizeSize(2)
                    if ~regionMask(newRow, newCol)
                        hueValue = imgHSV(newRow, newCol, 1);
                        saturationValue = imgHSV(newRow, newCol, 2);
                        brightnessValue = imgHSV(newRow, newCol, 3);
                        
                        % Add point to checklist regardless of criteria
                        checklist = [checklist; hueValue, saturationValue, brightnessValue];

                        % Check if the point fits the cell criteria
                        if hueValue >= hueLowThreshold && hueValue <= hueHighThreshold && ...
                           brightnessValue >= brightnessLowThreshold && brightnessValue <= brightnessHighThreshold && ...
                           saturationValue > saturationThreshold
                           avgHue = mean(checklist(end-9:end, 1));
                            
                            % Perform similarity check if average hue is below threshold
                            if avgHue < hueBackgroundThreshold
                                newHsvValue = [hueValue, saturationValue, brightnessValue];
                                avgSaturation = mean(hsvValues, 2);
                                hsvDistance = sqrt(sum((saturationValue - avgSaturation).^2));

                                if hsvDistance <= maxHsvDistance
                                    regionMask(newRow, newCol) = true;
                                    growingList = [growingList; newRow, newCol];
                                    hsvValues = [hsvValues; newHsvValue];
                                end
                            else
                                regionMask(newRow, newCol) = true;
                                growingList = [growingList; newRow, newCol];
                                hsvValues = [hsvValues; hueValue, saturationValue, brightnessValue];
                            end
                        end
                    end
                end
            end
        end
    end
    
    % Clean up the mask using morphological operations
    cellMask = imfill(regionMask, 'holes');
    cellMask = imopen(cellMask, strel('disk', 5));
    cellMask = imclose(cellMask, strel('disk', 5));
    
    % Create overlay image with segmentation
    imgOverlay = imgResized;
    imgOverlay(:, :, 1) = imgResized(:, :, 1) .* uint8(~cellMask) + uint8(cellMask) * 255;
    imgOverlay(:, :, 2) = imgResized(:, :, 2) .* uint8(~cellMask) + uint8(cellMask) * 0;
    imgOverlay(:, :, 3) = imgResized(:, :, 3) .* uint8(~cellMask) + uint8(cellMask) * 0;
    
    % Create a new figure for each cell type
    figure;
    % Plot results in a 1x3 subplot
    subplot(1, 3, 1); imshow(imgResized); title([cellType ' - Original']);
    subplot(1, 3, 2); imshow(cellMask); title([cellType ' - Cell Mask']);
    subplot(1, 3, 3); imshow(imgOverlay); title([cellType ' - Segmentation']);
    
    % Adjust figure size and layout
    set(gcf, 'Position', [100, 100, 1600, 400]);
end
