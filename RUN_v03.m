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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SEGMENTATION %%%%%%%%%%%%%

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
    maxHsvDistance = 0.4; % Maximum allowable HSV distance for new points
    
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
                           
                            % Check if the last 10 points fit the HSV criteria
                            if all(checklist(end-9:end, 1) >= hueLowThreshold) && ...
                               mean(checklist(end-9:end, 1)) >= hueBackgroundThreshold
                                % If all last 10 points have valid hue values, proceed as before
                                regionMask(newRow, newCol) = true;
                                growingList = [growingList; newRow, newCol];
                                hsvValues = [hsvValues; hueValue, saturationValue, brightnessValue];
                            else
                                % Calculate similarity based on saturation distance
                                avgSaturation = mean(hsvValues(:, 2));
                                hsvDistance = abs(saturationValue - avgSaturation);
    
                                % Perform similarity check based on HSV distance
                                if hsvDistance <= maxHsvDistance
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FEATURE EXTRACTION %%%%%%%%%%%%%

    % Define arrays to store extracted features
    numFeatures = 20;  % Number of features to extract (adjusted for hueHistogram)
    featureMatrix = zeros(length(cellTypes), numFeatures);  % Store features for each cell type
    
    % Loop through each cell type folder
    for i = 1:numCellTypes
        cellType = cellTypes{i};
        imageDir = fullfile(mainDir, cellType);
        imageFiles = dir(fullfile(imageDir, '*.jpg'));
        
        % Initialize feature vector for the current cell type
        features = zeros(length(imageFiles), numFeatures);  % Adjusted for multiple images per cell type
        
        % Process each image in the current cell type folder
        for j = 1:length(imageFiles)
            imgPath = fullfile(imageDir, imageFiles(j).name);
            img = imread(imgPath);
            
            % Resize image to the predefined size
            imgResized = imresize(img, resizeSize);
            imgHSV = rgb2hsv(imgResized);
            
            % Segment the image to obtain cell mask (already done in your code)
            % Replace or use 'cellMask' variable as per your implementation
            
            % Clean up the mask using morphological operations (already done in your code)
            cellMask = imfill(regionMask, 'holes');
            cellMask = imopen(cellMask, strel('disk', 5));
            cellMask = imclose(cellMask, strel('disk', 5));
            
            % Extract region properties from the cell mask
            props = regionprops(cellMask, 'Area', 'MajorAxisLength', 'MinorAxisLength', ...
                'EquivDiameter', 'Perimeter', 'Orientation', 'Solidity');
            
            % Calculate additional properties if needed
            area = props.Area;
            diameter = props.EquivDiameter;
            roundness = 4 * pi * area / props.Perimeter^2;  % Roundness measure
            
            % Calculate average hue within the cell mask
            hueValues = imgHSV(:,:,1);
            maskedHues = hueValues(cellMask);
            avgHue = mean(maskedHues);
            
            % Calculate histogram of hue values within the cell mask
            hueHistogram = histcounts(maskedHues, linspace(0, 1, numFeatures));  % Adjusted to numFeatures bins
            
            % Store features in the feature vector
            features(j, 1) = area;
            features(j, 2) = diameter;
            features(j, 3) = roundness;
            features(j, 4) = avgHue;
            features(j, 5:numFeatures) = hueHistogram;  % Assign all elements of hueHistogram
            
            % Store features for this image
            featureMatrix(i, :) = mean(features, 1);  % Average features across images for this cell type
        end
        
        % Display features extracted for the current cell type
        disp(['Features for ' cellType ':']);
        disp(featureMatrix(i, :));
    end
    
    % Display the entire feature matrix for review
    disp('Feature Matrix:');
    disp(featureMatrix);
end
