% RUN.m - Main script to classify different cell types on microscopy images

% Clear workspace and close all figures
clear; close all; clc;

% Set random seed for reproducibility
rng(42);

% Define the main directory containing images of different cell types
mainDir = 'C:\Users\Daniel\Documents\FHKärnten\Medical Engineering\2. Semester\Applied Medical Image Analysis\Project_Part2\bloodcells_4classes_500_test';
maskFolder = 'C:\Users\Daniel\Documents\FHKärnten\Medical Engineering\2. Semester\Applied Medical Image Analysis\Project_Part2\bloodcells_4classes_500_test\mask';
cellTypes = {'basophil', 'eosinophil', 'erythroblast', 'ig'};
numCellTypes = length(cellTypes);

% Define image size to resize to
resizeSize = [256, 256];

% Initialize an empty cell array to store features with headers
headers = {'Area', 'Diameter', 'Circularity', 'Perimeter', 'MajorAxisLength', 'MinorAxisLength', ...
           'Eccentricity', 'Curvature', 'Roughness', 'Avg Hue', 'Avg Saturation', 'Avg Brightness', ...
           'Hue 0.6+', 'Hue 0.7+', 'Hue 0.8+', 'Hue 0.9+', ...
           'Number Superpixels', 'Number Full Edges', 'Longest Inner Edge'};
featureMatrix = headers;
featureMatrixArray = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SEGMENTATION %%%%%%%%%%%%%

% Loop through each cell type folder
for i = 1:numCellTypes
    featureMatrixArray = [];
    cellType = cellTypes{i};
    imageDir = fullfile(mainDir, cellType);
    imageFiles = dir(fullfile(imageDir, '*.jpg'));
    features = [];
    
    % Process each image in the current cell type folder
    for j = 1:length(imageFiles)

        imgPath = fullfile(imageDir, imageFiles(j).name);
        img = imread(imgPath);
        
        % Resize image to the predefined size
        imgResized = imresize(img, resizeSize);
        imgHSV = rgb2hsv(imgResized);
        numPixels = resizeSize(1)*resizeSize(2);
        
        % Convert image to grayscale for visualization (optional)
        imgGray = rgb2gray(imgResized);
        
        % Define the seed point (center of the image)
        seedPoint = [round(resizeSize(1)/2), round(resizeSize(2)/2)];
        
        % Region growing parameters (example values, adjust as needed)
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
        
        % Create overlay image with segmentation (optional)
        imgOverlay = imgResized;
        imgOverlay(:, :, 1) = imgResized(:, :, 1) .* uint8(~cellMask) + uint8(cellMask) * 255;
        imgOverlay(:, :, 2) = imgResized(:, :, 2) .* uint8(~cellMask) + uint8(cellMask) * 0;
        imgOverlay(:, :, 3) = imgResized(:, :, 3) .* uint8(~cellMask) + uint8(cellMask) * 0;
        
        % Display the segmented image (optional)
        figure;
        subplot(1, 3, 1); imshow(imgResized); title([cellType ' - Original']);
        subplot(1, 3, 2); imshow(cellMask); title([cellType ' - Cell Mask']);
        subplot(1, 3, 3); imshow(imgOverlay); title([cellType ' - Segmentation']);
        set(gcf, 'Position', [100, 100, 1600, 400]);
        
       
%%%%%%%%%%%%%%%% FEATURE EXTRACTION %%%%%%%%%%%%%%%%%%%%%

        % Perform superpixel segmentation as feature
        numSuperpixels = 8;  % Example: Number of superpixels
        SuperPixelInput = double(imgGray) .* cellMask;
        [superpixelLabels, numSuperpixels] = superpixels(SuperPixelInput, numSuperpixels, "Compactness",5, "Method",'slic');
        % Find unique superpixel labels within the segmented region
        uniqueSuperpixelLabels = unique(superpixelLabels(cellMask));
        
        % Count the number of unique superpixels within the segmented region
        numUniqueSuperpixels = numel(uniqueSuperpixelLabels);
        % Create RGB overlay image
        superpixelOverlay = label2rgb(superpixelLabels);
        figure;
        imshow(imgGray);  % Display hue channel
        hold on;
        h = imshow(superpixelOverlay);
        set(h, 'AlphaData', 0.5);  % Set transparency to 50%
        title('Superpixel Overlay on Gray Channel');
        
         % Extract region properties from the superpixels mask
        superpixelprops = regionprops(numSuperpixels, 'Area', 'Perimeter', 'Circularity','PixelIdxList');
        for k = 1:numUniqueSuperpixels
            boundary = bwboundaries(superpixelLabels == k);
            superarea(k) = superpixelprops(k).Area; 
            superperimeter(k) = superpixelprops(k).Perimeter;
            supercirc(k) = superpixelprops(k).Circularity;
        end



        % Perform Canny edge detection as feature
        blurredImg = imgaussfilt(double(imgGray) .* cellMask, 1.5);
        edges = edge(blurredImg .* cellMask, 'Sobel');
        
        % Display the Canny edges
        figure;
        imshow(edges);
        title('Canny Edges of imgGray .* cellMask');
        
        % Count number of full (strong) edges
        numEdges = sum(edges(:));

        disp(['Number of full edges (strong edges): ' num2str(numEdges)]);
        % Find connected components (continuous edges) in the binary edge image
        cc = bwconncomp(edges);
        
        % Measure the lengths of the connected components (edges)
        edgeLengths = cellfun(@length, cc.PixelIdxList);
        numLongEdges = edgeLengths > 20;
        % Find the (second/third) longest edge
        [sortedLengths, sortedIdx] = sort(edgeLengths, 'descend');
        



        % Extract region properties from the cell mask
        props = regionprops(cellMask, 'Area', 'Circularity', 'MajorAxisLength', 'MinorAxisLength', ...
            'EquivDiameter', 'Perimeter', 'Orientation', 'Solidity', 'MajorAxisLength', 'MinorAxisLength', 'Eccentricity');
        
        % Calculate additional properties if needed
        areapixel = props.Area;
        area = areapixel / numPixels;
        diameter = props.EquivDiameter;
        circularity = props.Circularity;
        perimeter = props.Perimeter;
        majoraxis = props.MajorAxisLength;
        minoraxis = props.MinorAxisLength;
        eccentricity = props.Eccentricity;
        curvature = (perimeter^2) / (4 * pi * areapixel);  % Curvature measure
        roughness = perimeter / areapixel;  % Roughness measure
        
        % Calculate average properties of imgHSV within the cell mask
        maskedHue = imgHSV(:,:,1) .* cellMask;
        maskedSaturation = imgHSV(:,:,2) .* cellMask;
        maskedBrightness = imgHSV(:,:,3) .* cellMask;
        
        % Calculate average hue, saturation, and brightness within the cell mask
        avgHue = mean(maskedHue(maskedHue > 0), 'all', 'omitnan');
        avgSaturation = mean(maskedSaturation(maskedSaturation > 0), 'all', 'omitnan');
        avgBrightness = mean(maskedBrightness(maskedBrightness > 0), 'all', 'omitnan');
        
        % Additional features: number of pixels within specific hue ranges
        % Define hue ranges
        hueRanges = [0.6, 0.7, 0.8, 0.9, 1.0];
        numPixelsInRange = zeros(1, length(hueRanges)-1);
        
        for k = 1:length(hueRanges)-1
            lowerHue = hueRanges(k);
            upperHue = hueRanges(k+1);
            numPixelsInRange(k) = sum(maskedHue >= lowerHue & maskedHue < upperHue, 'all', 'omitnan');
            areaPixelsInRange(k) = numPixelsInRange(k)/areapixel;
        end
        
        % Store features in the feature vector
        features(j, 1) = area;
        features(j, 2) = diameter;
        features(j, 3) = circularity;
        features(j, 4) = perimeter;
        features(j, 5) = majoraxis;
        features(j, 6) = minoraxis;
        features(j, 7) = eccentricity;
        features(j, 8) = curvature;
        features(j, 9) = roughness;
        features(j, 10) = avgHue;
        features(j, 11) = avgSaturation;
        features(j, 12) = avgBrightness;
        features(j, 13:16) = areaPixelsInRange;
        features(j, 17) = numUniqueSuperpixels;
        features(j, 18) = size(numLongEdges,1);
        features(j, 19) = sortedLengths(2) / perimeter;
        
        % Display features extracted for the current image (optional)
        disp(['Features for ' cellType ' - Image ' num2str(j) ':']);
        disp(features(j, :));
    end
    featureMatrix = [featureMatrix; num2cell(features)];
    featureMatrixArray = [featureMatrixArray; features];
    % Save the feature matrix for the current cell type
    saveFileName = fullfile(mainDir, [cellType '_features.mat']);
    save(saveFileName, 'featureMatrix');
    % Save the feature matrix as a CSV file
    saveFileNameCSV = fullfile(mainDir, [cellType '_features.csv']);
    csvwrite(saveFileNameCSV, featureMatrixArray);



end


















    
% end
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FEATURE EXTRACTION %%%%%%%%%%%%%
%  % Define arrays to store extracted features
% numFeatures = 16;  % Number of features to extract
% featureMatrix = zeros(numCellTypes, numFeatures);  % Store features for each cell type
% 
% % Loop through each cell type folder
% for i = 1:numCellTypes
%     cellType = cellTypes{i};
%     imageDir = fullfile(mainDir, cellType);
%     imageFiles = dir(fullfile(imageDir, '*.jpg'));
% 
%     % Initialize feature vector for the current cell type
%     features = zeros(length(imageFiles), numFeatures);
% 
%     % Process each image in the current cell type folder
%     for j = 1:length(imageFiles)
%         imgPath = fullfile(imageDir, imageFiles(j).name);
%         img = imread(imgPath);
% 
%         % Resize image to the predefined size
%         imgResized = imresize(img, resizeSize);
%         imgHSV = rgb2hsv(imgResized);
%         numPixels = resizeSize(1)*resizeSize(2);
% 
%         % Segment the image to obtain cell mask (already done in your segmentation step)
%         % Replace or use 'cellMask' variable as per your implementation
% 
%         % Assuming 'cellMask' is already obtained from segmentation
%         % Clean up the mask using morphological operations (already done in your segmentation step)
%         % For illustration, using some example values
%         cellMask = logical(randi([0 1], resizeSize));  % Replace with your actual cell mask
% 
%         % Extract region properties from the cell mask
%         props = regionprops(cellMask, 'Area', 'Circularity', 'MajorAxisLength', 'MinorAxisLength', ...
%             'EquivDiameter', 'Perimeter', 'Orientation', 'Solidity', 'MajorAxisLength', 'MinorAxisLength', 'Eccentricity');
% 
%         % Calculate additional properties if needed
%         area = props.Area;
%         area = area / numPixels;
%         diameter = props.EquivDiameter;
%         circularity = props.Circularity;
%         perimeter = props.Perimeter;
%         majoraxis = props.MajorAxisLength;
%         minoraxis = props.MinorAxisLength;
%         eccentricity = props.Eccentricity;
%         curvature = (perimeter^2) / (4 * pi * area);  % Curvature measure
%         roughness = perimeter / area;  % Roughness measure
% 
%         % Calculate average properties of imgHSV within the cell mask
%         maskedHue = imgHSV(:,:,1) .* cellMask;
%         maskedSaturation = imgHSV(:,:,2) .* cellMask;
%         maskedBrightness = imgHSV(:,:,3) .* cellMask;
% 
%         % Calculate average hue, saturation, and brightness within the cell mask
%         avgHue = mean(maskedHue(maskedHue > 0), 'all', 'omitnan');
%         avgSaturation = mean(maskedSaturation(maskedSaturation > 0), 'all', 'omitnan');
%         avgBrightness = mean(maskedBrightness(maskedBrightness > 0), 'all', 'omitnan');
% 
%         % Additional features: number of pixels within specific hue ranges
%         % Define hue ranges
%         hueRanges = [0.6, 0.7, 0.8, 0.9, 1.0];
%         numPixelsInRange = zeros(1, length(hueRanges)-1);
% 
%         for k = 1:length(hueRanges)-1
%             lowerHue = hueRanges(k);
%             upperHue = hueRanges(k+1);
%             numPixelsInRange(k) = sum(maskedHue >= lowerHue & maskedHue < upperHue, 'all', 'omitnan');
%             areaPixelsInRange(k) = numPixelsInRange(k)/area;
%         end
% 
%         % Store features in the feature vector
%         features(j, 1) = area;
%         features(j, 2) = diameter;
%         features(j, 3) = circularity;
%         features(j, 4) = perimeter;
%         features(j, 5) = majoraxis;
%         features(j, 6) = minoraxis;
%         features(j, 7) = eccentricity;
%         features(j, 8) = curvature;
%         features(j, 9) = roughness;
%         features(j, 10) = avgHue;
%         features(j, 11) = avgSaturation;
%         features(j, 12) = avgBrightness;
%         features(j, 13:16) = areaPixelsInRange;
% 
%         % Display features extracted for the current image (optional)
%         disp(['Features for ' cellType ' - Image ' num2str(j) ':']);
%         disp(features(j, :));
%     end
% 
%     % Store features for the current cell type in the feature matrix
%     featureMatrix(i, :) = mean(features, 1);  % Using mean of features across images
%     % Save the feature matrix for the current cell type
%     saveFileName = fullfile(mainDir, [cellType '_features.mat']);
%     save(saveFileName, 'featureMatrix');
%     % Save the feature matrix as a CSV file
%     saveFileNameCSV = fullfile(mainDir, [cellType '_features.csv']);
%     csvwrite(saveFileNameCSV, features);
% 
% end



% Display confirmation
disp(['Feature matrix for ' cellType ' saved as ' saveFileName]);

% Display the final feature matrix for all cell types (if needed)
disp('Final Feature Matrix:');
disp(featureMatrix);