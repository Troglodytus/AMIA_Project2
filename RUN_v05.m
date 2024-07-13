% RUN.m - Main script to classify different cell types on microscopy images

% Clear workspace and close all figures
clear; close all; clc;

% Set random seed for reproducibility
rng(42);

% Define the main directory containing images of different cell types
mainDir = 'C:\Users\Daniel\Documents\FHKärnten\Medical Engineering\2. Semester\Applied Medical Image Analysis\Project_Part2\bloodcells_4classes_20';
maskFolder = 'C:\Users\Daniel\Documents\FHKärnten\Medical Engineering\2. Semester\Applied Medical Image Analysis\Project_Part2\bloodcells_4classes_20\mask';
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
allFeatures = [];
allLabels = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SEGMENTATION %%%%%%%%%%%%%

% Loop through each cell type folder
for i = 1:numCellTypes
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
        
        % Append features and labels for the current image to the overall arrays
        allFeatures = [allFeatures; features(j, :)];
        allLabels = [allLabels; {cellType}];
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RANDOM FOREST TRAINING %%%%%%%%%%%%%

% Convert cell array of labels to a categorical array
allLabels = categorical(allLabels);

% Train a Random Forest classifier using TreeBagger
nTrees = 100;  % Number of trees in the forest
inBagFraction = 0.7;  % Fraction of data used for training each tree

% Create the TreeBagger object
rfModel = TreeBagger(nTrees, allFeatures, allLabels, ...
    'OOBPredictorImportance', 'on', ...
    'InBagFraction', inBagFraction);

% Display out-of-bag error
oobError = oobError(rfModel);
figure;
plot(oobError);
xlabel('Number of Grown Trees');
ylabel('Out-of-Bag Classification Error');
title('Out-of-Bag Error vs. Number of Grown Trees');

% Display predictor importance
predictorImportance = rfModel.OOBPermutedPredictorDeltaError;
figure;
bar(predictorImportance);
xlabel('Predictor Index');
ylabel('Predictor Importance');
title('Predictor Importance Estimates');

% Save the Random Forest model
save(fullfile(mainDir, 'rfModel.mat'), 'rfModel');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% APPLY CLASSIFIER TO TEST DATA %%%%%%%%%%%%%

% Assuming test data is structured similarly to training data
% and saved in a separate folder named 'test'

testDir = 'C:\Users\Daniel\Documents\FHKärnten\Medical Engineering\2. Semester\Applied Medical Image Analysis\Project_Part2\bloodcells_4classes_500_test\';
testCellTypes = {'basophil', 'eosinophil', 'erythroblast', 'ig'};
numTestCellTypes = length(testCellTypes);

testFeatures = [];
testLabels = [];

for i = 1:numTestCellTypes
    cellType = testCellTypes{i};
    testImageDir = fullfile(testDir, cellType);
    testImageFiles = dir(fullfile(testImageDir, '*.jpg'));
    
    for j = 1:length(testImageFiles)
        imgPath = fullfile(testImageDir, testImageFiles(j).name);
        img = imread(imgPath);
        
        % Resize image to the predefined size
        imgResized = imresize(img, resizeSize);
        imgHSV = rgb2hsv(imgResized);
        numPixels = resizeSize(1)*resizeSize(2);
        
        % Perform segmentation and feature extraction as done for training data
        % Initialize variables for feature extraction
        features = zeros(1, size(allFeatures, 2)); % Ensure the correct size
        
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
        
        % Perform superpixel segmentation as feature
        numSuperpixels = 8;  % Example: Number of superpixels
        SuperPixelInput = double(rgb2gray(imgResized)) .* cellMask;
        [superpixelLabels, numSuperpixels] = superpixels(SuperPixelInput, numSuperpixels, "Compactness",5, "Method",'slic');
        % Find unique superpixel labels within the segmented region
        uniqueSuperpixelLabels = unique(superpixelLabels(cellMask));
        
        % Count the number of unique superpixels within the segmented region
        numUniqueSuperpixels = numel(uniqueSuperpixelLabels);
        
        % Perform Canny edge detection as feature
        blurredImg = imgaussfilt(double(rgb2gray(imgResized)) .* cellMask, 1.5);
        edges = edge(blurredImg .* cellMask, 'Sobel');
        
        % Count number of full (strong) edges
        numEdges = sum(edges(:));

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
        avgHue = mean(maskedHue(maskedHue > 0));
        avgSaturation = mean(maskedSaturation(maskedSaturation > 0));
        avgBrightness = mean(maskedBrightness(maskedBrightness > 0));
        
        % Set features
        features(1) = area;
        features(2) = diameter;
        features(3) = circularity;
        features(4) = perimeter;
        features(5) = majoraxis;
        features(6) = minoraxis;
        features(7) = eccentricity;
        features(8) = curvature;
        features(9) = roughness;
        features(10) = avgHue;
        features(11) = avgSaturation;
        features(12) = avgBrightness;
        features(13) = numUniqueSuperpixels;
        features(14) = numEdges;
        
        % Store features and labels for the current test image
        testFeatures = [testFeatures; features];
        testLabels = [testLabels; {cellType}];
    end
end

% Convert cell array of test labels to a categorical array
testLabels = categorical(testLabels);

% Check if the number of columns in testFeatures matches the training features
if size(testFeatures, 2) ~= size(allFeatures, 2)
    error('The number of columns in testFeatures does not match the training features.');
end

% Predict using the trained Random Forest model
[predictedLabels, scores] = predict(rfModel, testFeatures);
% Convert predicted labels to categorical array to match testLabels
predictedLabels = categorical(predictedLabels);

% Calculate accuracy
accuracy = sum(predictedLabels == testLabels) / length(testLabels);
disp(['Test Accuracy: ' num2str(accuracy * 100) '%']);

% Display confusion matrix
confMat = confusionmat(testLabels, predictedLabels);
figure;
confusionchart(confMat, testCellTypes);
title('Confusion Matrix for Test Data');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DISPLAY TREE STRUCTURE %%%%%%%%%%%%%

% Display structure of three trees
figure;
view(rfModel.Trees{1}, 'Mode', 'graph');
title('Structure of Tree 1');

figure;
view(rfModel.Trees{2}, 'Mode', 'graph');
title('Structure of Tree 2');

figure;
view(rfModel.Trees{3}, 'Mode', 'graph');
title('Structure of Tree 3');

% Show the first of each category of images during training
if j == 1
    figure;
    subplot(2,2,1); imshow(imgResized); title('Original Resized Image');
    subplot(2,2,2); imshow(cellMask); title('Cell Mask');
    subplot(2,2,3); imshow(label2rgb(superpixelLabels)); title('Superpixels');
    subplot(2,2,4); imshow(edges); title('Edges');
end

