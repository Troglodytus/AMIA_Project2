% RUN.m - Main script to classify different cell types on microscopy images

% Clear workspace and close all figures
clear; close all; clc;

% Set random seed for reproducibility
rng(42);

% Define the main directory containing images of different cell types
mainDir = 'C:\Users\Daniel\Documents\FHKärnten\Medical Engineering\2. Semester\Applied Medical Image Analysis\Project_Part2\bloodcells_4classes_500';
maskFolder = 'C:\Users\Daniel\Documents\FHKärnten\Medical Engineering\2. Semester\Applied Medical Image Analysis\Project_Part2\bloodcells_4classes_500\mask';
cellTypes = {'basophil', 'eosinophil', 'erythroblast', 'ig'};
numCellTypes = length(cellTypes);

% Define image size to resize to
resizeSize = [256, 256];

% Initialize an empty cell array to store features with headers
headers = {'Area', 'Diameter', 'Circularity', 'Perimeter', 'MajorAxisLength', 'MinorAxisLength', ...
           'Eccentricity', 'Curvature', 'Roughness', 'Avg Hue', 'Avg Saturation', 'Avg Brightness', 'Kurtosis Intensity', 'Cell Entropy', 'Hist Frequency Median','Cell Gradient Magnitude', ...
           'Hue 0.6+', 'Hue 0.7+', 'Hue 0.8+', 'Hue 0.9+', ...
           'Number Superpixels', 'Number Full Edges', 'Longest Inner Edge'};
featureMatrix = headers;
featureMatrixArray = [];
allFeatures = [];
allLabels = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SEGMENTATION %%%%%%%%%%%%%

% % Loop through each cell type folder
% for i = 1:numCellTypes
%     cellType = cellTypes{i};
%     imageDir = fullfile(mainDir, cellType);
%     imageFiles = dir(fullfile(imageDir, '*.jpg'));
%     features = [];
% 
%     % Process each image in the current cell type folder
%     for j = 1:length(imageFiles)
% 
%         imgPath = fullfile(imageDir, imageFiles(j).name);
%         img = imread(imgPath);
% 
%         % Resize image to the predefined size
%         imgResized = imresize(img, resizeSize);
%         imgHSV = rgb2hsv(imgResized);
%         numPixels = resizeSize(1)*resizeSize(2);
% 
%         % Convert image to grayscale for visualization (optional)
%         imgGray = rgb2gray(imgResized);
% 
%         % Define the seed point (center of the image)
%         seedPoint = [round(resizeSize(1)/2), round(resizeSize(2)/2)];
% 
%         % Region growing parameters (example values, adjust as needed)
%         hueLowThreshold = 0.05;
%         hueHighThreshold = 0.35;
%         brightnessLowThreshold = 0.24;
%         brightnessHighThreshold = 0.84;
%         saturationThreshold = 0.1;
%         hueBackgroundThreshold = 0.5;
%         maxHsvDistance = 0.32; % Maximum allowable HSV distance for new points
% 
%         % Initialize region growing
%         regionMask = false(resizeSize);
%         regionMask(seedPoint(1), seedPoint(2)) = true;
%         growingList = seedPoint;
%         hsvValues = squeeze(imgHSV(seedPoint(1), seedPoint(2), :))';
%         checklist = hsvValues; % Initialize checklist with the seed point HSV values
% 
%         % Perform region growing for initially 10 steps
%         for step = 1:10
%             if isempty(growingList)
%                 break;
%             end
% 
%             currentPoint = growingList(1, :);
%             growingList(1, :) = [];
% 
%             for k = -1:1
%                 for l = -1:1
%                     newRow = currentPoint(1) + k;
%                     newCol = currentPoint(2) + l;
% 
%                     if newRow > 0 && newRow <= resizeSize(1) && newCol > 0 && newCol <= resizeSize(2)
%                         if ~regionMask(newRow, newCol)
%                             hueValue = imgHSV(newRow, newCol, 1);
%                             saturationValue = imgHSV(newRow, newCol, 2);
%                             brightnessValue = imgHSV(newRow, newCol, 3);
% 
%                             % Add point to checklist regardless of criteria
%                             checklist = [checklist; hueValue, saturationValue, brightnessValue];
% 
%                             % Check if the point fits the cell criteria
%                             if hueValue <= hueLowThreshold || hueValue >= hueHighThreshold && ...
%                                     brightnessValue >= brightnessLowThreshold && brightnessValue <= brightnessHighThreshold && ...
%                                     saturationValue > saturationThreshold
%                                 regionMask(newRow, newCol) = true;
%                                 growingList = [growingList; newRow, newCol];
%                                 hsvValues = [hsvValues; hueValue, saturationValue, brightnessValue];
%                             end
%                         end
%                     end
%                 end
%             end
%         end
% 
%         % Perform the main region growing with additional similarity check
%         while ~isempty(growingList)
%             currentPoint = growingList(1, :);
%             growingList(1, :) = [];
% 
%             for k = -1:1
%                 for l = -1:1
%                     newRow = currentPoint(1) + k;
%                     newCol = currentPoint(2) + l;
% 
%                     if newRow > 0 && newRow <= resizeSize(1) && newCol > 0 && newCol <= resizeSize(2)
%                         if ~regionMask(newRow, newCol)
%                             hueValue = imgHSV(newRow, newCol, 1);
%                             saturationValue = imgHSV(newRow, newCol, 2);
%                             brightnessValue = imgHSV(newRow, newCol, 3);
% 
%                             % Add point to checklist regardless of criteria
%                             checklist = [checklist; hueValue, saturationValue, brightnessValue];
% 
%                             % Check if the point fits the cell criteria
%                             if hueValue <= hueLowThreshold || hueValue >= hueHighThreshold && ...
%                                     brightnessValue >= brightnessLowThreshold && brightnessValue <= brightnessHighThreshold && ...
%                                     saturationValue > saturationThreshold
% 
%                                 % Check if the last 10 points fit the HSV criteria
%                                 if all(checklist(end-9:end, 1) >= hueHighThreshold) && ...
%                                         mean(checklist(end-9:end, 1)) >= hueBackgroundThreshold
%                                     regionMask(newRow, newCol) = true;
%                                     growingList = [growingList; newRow, newCol];
%                                     hsvValues = [hsvValues; hueValue, saturationValue, brightnessValue];
%                                 else
%                                     % Calculate similarity based on saturation distance
%                                     avgSaturation = mean(hsvValues(:, 2));
%                                     hsvDistance = abs(saturationValue - avgSaturation);
% 
%                                     % Perform similarity check based on HSV distance
%                                     if hsvDistance <= maxHsvDistance
%                                         regionMask(newRow, newCol) = true;
%                                         growingList = [growingList; newRow, newCol];
%                                         hsvValues = [hsvValues; hueValue, saturationValue, brightnessValue];
%                                     end
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%         end
% 
%         % Clean up the mask using morphological operations
%         cellMask = imfill(regionMask, 'holes');
%         cellMask = imopen(cellMask, strel('disk', 5));
%         cellMask = imclose(cellMask, strel('disk', 5));
% 
%         % Create overlay image with segmentation
%         imgOverlay = imgResized;
%         imgOverlay(:, :, 1) = imgResized(:, :, 1) .* uint8(~cellMask) + uint8(cellMask) * 255;
%         imgOverlay(:, :, 2) = imgResized(:, :, 2) .* uint8(~cellMask) + uint8(cellMask) * 0;
%         imgOverlay(:, :, 3) = imgResized(:, :, 3) .* uint8(~cellMask) + uint8(cellMask) * 0;
% 
%         % Display the segmented image (optional)
%         % figure;
%         % subplot(1, 3, 1); imshow(imgResized); title([cellType ' - Original']);
%         % subplot(1, 3, 2); imshow(cellMask); title([cellType ' - Cell Mask']);
%         % subplot(1, 3, 3); imshow(imgOverlay); title([cellType ' - Segmentation']);
%         % set(gcf, 'Position', [100, 100, 1600, 400]);
% 
% 
% %%%%%%%%%%%%%%%% FEATURE EXTRACTION %%%%%%%%%%%%%%%%%%%%%
% 
%         % Perform superpixel segmentation as feature
%         numSuperpixels = 8;  % Example: Number of superpixels
%         SuperPixelInput = double(imgGray) .* cellMask;
%         [superpixelLabels, numSuperpixels] = superpixels(SuperPixelInput, numSuperpixels, "Compactness",5, "Method",'slic');
%         % Find unique superpixel labels within the segmented region
%         uniqueSuperpixelLabels = unique(superpixelLabels(cellMask));
% 
%         % Count the number of unique superpixels within the segmented region
%         numUniqueSuperpixels = numel(uniqueSuperpixelLabels);
%         % Create RGB overlay image
%         superpixelOverlay = label2rgb(superpixelLabels);
%         % Display the Superpixels (optional)
%         % figure;
%         % imshow(imgGray);  % Display hue channel
%         % hold on;
%         % h = imshow(superpixelOverlay);
%         % set(h, 'AlphaData', 0.5);  % Set transparency to 50%
%         % title('Superpixel Overlay on Gray Channel');
% 
%          % Extract region properties from the superpixels mask
%         superpixelprops = regionprops(numSuperpixels, 'Area', 'Perimeter', 'Circularity','PixelIdxList');
%         for k = 1:numUniqueSuperpixels
%             boundary = bwboundaries(superpixelLabels == k);
%             superarea(k) = superpixelprops(k).Area; 
%             superperimeter(k) = superpixelprops(k).Perimeter;
%             supercirc(k) = superpixelprops(k).Circularity;
%         end
% 
% 
% 
%         % Perform Canny edge detection as feature
%         blurredImg = imgaussfilt(double(imgGray) .* cellMask, 1.5);
%         edges = edge(blurredImg .* cellMask, 'Sobel');
% 
%         % Display the Canny edge image (optional)
%         % figure;
%         % imshow(edges);
%         % title('Canny Edges of imgGray .* cellMask');
% 
%         % Count number of full (strong) edges
%         numEdges = sum(edges(:));
%         %disp(['Number of full edges (strong edges): ' num2str(numEdges)]);
% 
%         % Find connected components (continuous edges) in the binary edge image
%         cc = bwconncomp(edges);
% 
%         % Measure the lengths of the connected components (edges)
%         edgeLengths = cellfun(@length, cc.PixelIdxList);
%         numLongEdges = edgeLengths > 20;
%         % Find the (second/third) longest edge
%         [sortedLengths, sortedIdx] = sort(edgeLengths, 'descend');
% 
% 
% 
% 
%         % Extract region properties from the cell mask
%         props = regionprops(cellMask, 'Area', 'Circularity', 'MajorAxisLength', 'MinorAxisLength', ...
%             'EquivDiameter', 'Perimeter', 'Orientation', 'Solidity', 'MajorAxisLength', 'MinorAxisLength', 'Eccentricity');
% 
%         % Check if all necessary properties are available
%         if isempty(props) || ~isfield(props, 'Area') || ~isfield(props, 'Circularity') || ...
%                 ~isfield(props, 'MajorAxisLength') || ~isfield(props, 'MinorAxisLength') || ...
%                 ~isfield(props, 'EquivDiameter') || ~isfield(props, 'Perimeter') || ...
%                 ~isfield(props, 'Orientation') || ~isfield(props, 'Solidity') || ...
%                 ~isfield(props, 'Eccentricity')
%             % Skip this image if any required properties are missing
%             disp(['Skipping image ' num2str(j) ' in folder ' cellType ': Missing required properties.']);
%             continue;  % Skip to the next image
%         end
% 
% 
%         % Calculate additional properties if needed
%         areapixel = props.Area;
%         area = areapixel / numPixels;
%         diameter = props.EquivDiameter;
%         circularity = props.Circularity;
%         perimeter = props.Perimeter;
%         majoraxis = props.MajorAxisLength;
%         minoraxis = props.MinorAxisLength;
%         eccentricity = props.Eccentricity;
%         curvature = (perimeter^2) / (4 * pi * areapixel);  % Curvature measure
%         roughness = perimeter / areapixel;  % Roughness measure
% 
%         % Calculate average properties of imgHSV within the cell mask
%         maskedHue = imgHSV(:,:,1) .* cellMask;
%         maskedSaturation = imgHSV(:,:,2) .* cellMask;
%         maskedBrightness = imgHSV(:,:,3) .* cellMask;
%         maskedIntensity = SuperPixelInput;
%         KurtosisIntensity = kurtosis(maskedIntensity(cellMask));
% 
%         %Entropy
%         maskedIntensityValues = maskedIntensity(cellMask & maskedIntensity > 0);
% 
%         %histogram of masked area
%         [counts, binLocations] = imhist(uint8(maskedIntensityValues * 255 / max(maskedIntensityValues)));
% 
%         % Normalize the histogram counts to get probabilities
%         p = counts / sum(counts);
% 
%         % Remove zero probabilities
%         p = p(p > 0);
% 
%         % Calculate entropy
%         cellEntropyIntensity = -sum(p .* log2(p));
% 
%         % Calculate the cumulative distribution function (CDF)
%         cdf = cumsum(counts) / sum(counts);
% 
%         % Find the intensity value at the median frequency
%         medianFrequencyValue = binLocations(find(cdf >= 0.5, 1));
% 
%         % Compute the gradient magnitude of the masked intensity
%         [~, gmag] = imgradient(maskedIntensity);
%         maskedGradient = gmag .* cellMask;
%         meanGradientMagnitude = mean(maskedGradient(maskedGradient > 0));
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
%             areaPixelsInRange(k) = numPixelsInRange(k)/areapixel;
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
%         features(j, 13) = KurtosisIntensity;
%         features(j, 14) = cellEntropyIntensity;
%         features(j, 15) = medianFrequencyValue;
%         features(j, 16) = meanGradientMagnitude;
%         features(j, 17:20) = areaPixelsInRange;
%         features(j, 21) = numUniqueSuperpixels;
%         features(j, 22) = size(numLongEdges,1);
%         if numLongEdges > 1
%             features(j, 23) = sortedLengths(2) / perimeter;
%         else
%             features(j, 23) = sortedLengths(1) / perimeter;
%         end
% 
% 
%         % Display features extracted for the current image (optional)
%         %disp(['Features for ' cellType ' - Image ' num2str(j) ':']);
%         %disp(features(j, :));
% 
%         % Append features and labels for the current image to the overall arrays
%         allFeatures = [allFeatures; features(j, :)];
%         allLabels = [allLabels; {cellType}];
%     end
% 
%     featureMatrix = [featureMatrix; num2cell(features)];
%     featureMatrixArray = [featureMatrixArray; features];
%     featureMatrixArray = allFeatures;
% 
%     % Save the feature matrix for the current cell type
%     saveFileName = fullfile(mainDir, [cellType '_features.mat']);
%     save(saveFileName, 'featureMatrix');
% 
%     %Save Labels
%     save(fullfile(mainDir, 'allLabels.mat'), 'allLabels');
%     %Save Feature Matrix
%     save(fullfile(mainDir, 'featureMatrixArray.mat'), 'featureMatrixArray');
% 
%     % Save the feature matrix as a CSV file
%     saveFileNameCSV = fullfile(mainDir, [cellType '_features.csv']);
%     csvwrite(saveFileNameCSV, featureMatrixArray);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%% LOAD FEATURES & DETECT OUTLIERS %%%%%%%%%%%%%%%%%%%%%%%%%%

% Load the data
load(fullfile(mainDir, 'allLabels.mat'), 'allLabels');
load(fullfile(mainDir, 'featureMatrixArray.mat'), 'featureMatrixArray');

% Check the data type of allLabels
if iscell(allLabels)
    allLabels = categorical(allLabels);
end

% Define a threshold for identifying outliers (e.g., Z-score > 3)
zScoreThreshold = 3;

% Initialize arrays to store filtered features and labels
filteredFeatures = [];
filteredLabels = categorical();

% Get unique cell types
cellTypes = unique(allLabels);

% Loop through each cell type
for i = 1:length(cellTypes)
    cellType = cellTypes(i);

    % Get indices of current cell type
    cellTypeIdx = (allLabels == cellType);

    % Extract features for the current cell type
    cellTypeFeatures = featureMatrixArray(cellTypeIdx, :);
    
    % Extract labels for the current cell type
    cellTypeLabels = allLabels(cellTypeIdx);

    % Calculate Z-scores for the features
    zScores = zscore(cellTypeFeatures);

    % Find rows without any outliers
    nonOutlierIdx = all(abs(zScores) <= zScoreThreshold, 2);

    % Append non-outlier features and labels to the filtered arrays
    filteredFeatures = [filteredFeatures; cellTypeFeatures(nonOutlierIdx, :)];
    filteredLabels = [filteredLabels; cellTypeLabels(nonOutlierIdx)];
end

% Save the filtered data
save(fullfile(mainDir, 'filtered_allLabels.mat'), 'filteredLabels');
save(fullfile(mainDir, 'filtered_featureMatrixArray.mat'), 'filteredFeatures');

% Save the filtered feature matrix as a CSV file
filteredFileNameCSV = fullfile(mainDir, 'filtered_featureMatrixArray.csv');
csvwrite(filteredFileNameCSV, filteredFeatures);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RANDOM FOREST TRAINING %%%%%%%%%%%%%

% Convert cell array of labels to a categorical array
filteredLabels = categorical(filteredLabels);

% Split data into training and test sets (80% training, 20% test)
cv = cvpartition(length(filteredLabels), 'HoldOut', 0.2);
idxTrain = cv.training;
idxTest = cv.test;

% Training data
trainFeatures = filteredFeatures(idxTrain, :);
trainLabels = filteredLabels(idxTrain);

% Test data
testFeatures = filteredFeatures(idxTest, :);
testLabels = filteredLabels(idxTest);

% Train a Random Forest classifier using TreeBagger
nTrees = 5000;  % Number of trees in the forest
inBagFraction = 0.75;  % Fraction of data used for training each tree

% Create the TreeBagger object
tic;
rfModel = TreeBagger(nTrees, trainFeatures, trainLabels, ...
    'OOBPredictorImportance', 'on', ...
    'InBagFraction', inBagFraction);

%Predict
predictedLabels = predict(rfModel, testFeatures);
classificationTimeFull = toc;
predictedLabels = categorical(predictedLabels);

% Compute accuracy
accuracy = sum(predictedLabels == testLabels) / numel(testLabels);
disp(['Accuracy on test set: ' num2str(accuracy * 100) '%']);
disp(['Classification time: ' num2str(classificationTimeFull) ' seconds']);

% Create confusion matrix
confMat = confusionmat(testLabels, predictedLabels);

% Display confusion matrix
figure;
heatmap(categories(filteredLabels), categories(filteredLabels), confMat);
xlabel('Predicted Labels');
ylabel('True Labels');
title('Confusion Matrix');

% Display out-of-bag error
oobErrors = oobError(rfModel);
figure;
plot(oobErrors);
xlabel('Number of Grown Trees');
ylabel('Out-of-Bag Classification Error');
title('Out-of-Bag Error vs. Number of Grown Trees');

% Display the out-of-bag classification error of each individual tree
numTrees = rfModel.NumTrees;
individualErrors = zeros(numTrees, 1);

for i = 1:numTrees
    % Get the out-of-bag predictions for the i-th tree
    [~, scores] = oobPredict(rfModel, 'Trees', i);
    oobSamples = find(rfModel.OOBIndices(:, i));  % Indices of OOB samples for the i-th tree
    oobPredictions = scores(oobSamples, :);       % Predictions for OOB samples
    [~, oobPredictedLabels] = max(oobPredictions, [], 2); % Convert scores to labels
    oobPredictedLabels = categorical(oobPredictedLabels, 1:length(categories(filteredLabels)), categories(filteredLabels));
    
    % Calculate OOB error for the i-th tree
    individualErrors(i) = sum(oobPredictedLabels ~= filteredLabels(oobSamples)) / length(oobSamples);
end

figure;
plot(1:numTrees, individualErrors, 'o-');
xlabel('Tree Index');
ylabel('Out-of-Bag Classification Error');
title('OOB Classification Error of Each Individual Tree');

% Display predictor importance
predictorImportance = rfModel.OOBPermutedPredictorDeltaError;
figure;
bar(predictorImportance);
xlabel('Predictor Index');
ylabel('Predictor Importance');
title('Predictor Importance Estimates');

% Display the structure of at least three trees
view(rfModel.Trees{1}, 'Mode', 'graph');
view(rfModel.Trees{2}, 'Mode', 'graph');
view(rfModel.Trees{3}, 'Mode', 'graph');

% Save the Random Forest model
save(fullfile(mainDir, 'rfModel.mat'), 'rfModel');


%%%%%%%%%%%%%%%%% PCA %%%%%%%%%%%%%%%%%%%%%%%%%

% Perform PCA on the feature matrix
[coeff, score, latent, tsquared, explained, mu] = pca(filteredFeatures);

explainedVarianceThreshold = 99.999; % Keep enough components to explain 99.5% of the variance
cumulativeVariance = cumsum(explained);
%numComponents = find(cumulativeVariance >= explainedVarianceThreshold, 1);

%Hardcode PCA Feature No.
numComponents = 18;

% Display number of components retained
disp(['Number of components retained: ', num2str(numComponents)]);

% Reduce the feature matrix to the selected number of principal components
reducedFeatures = score(:, 1:numComponents);


% Convert cell array of labels to a categorical array
filteredLabels = categorical(filteredLabels);

% Split data into training and test sets (80% training, 20% test)
cv = cvpartition(length(filteredLabels), 'HoldOut', 0.2);
idxTrain = cv.training;
idxTest = cv.test;

% Training data
trainFeatures = reducedFeatures(idxTrain, :);
trainLabels = filteredLabels(idxTrain);

% Test data
testFeatures = reducedFeatures(idxTest, :);
testLabels = filteredLabels(idxTest);

% Train a Random Forest classifier using TreeBagger
nTrees = 5000;  % Number of trees in the forest
inBagFraction = 0.75;  % Fraction of data used for training each tree

% Create the TreeBagger object after PCA
tic;
rfModelPCA = TreeBagger(nTrees, trainFeatures, trainLabels, ...
    'OOBPredictorImportance', 'on', ...
    'InBagFraction', inBagFraction);

% Predict on the test data
predictedLabels = predict(rfModelPCA, testFeatures);
classificationTimeFull = toc;
predictedLabels = categorical(predictedLabels);

% Compute accuracy
accuracy = sum(predictedLabels == testLabels) / numel(testLabels);
disp(['Accuracy on test set after PCA: ' num2str(accuracy * 100) '%']);
disp(['Classification time: ' num2str(classificationTimeFull) ' seconds']);

% Create confusion matrix
confMat = confusionmat(testLabels, predictedLabels);

% Display confusion matrix
figure;
heatmap(categories(filteredLabels), categories(filteredLabels), confMat);
xlabel('Predicted Labels');
ylabel('True Labels');
title('Confusion Matrix after PCA');

% Display out-of-bag error
oobErrors = oobError(rfModelPCA);
figure;
plot(oobErrors);
xlabel('Number of Grown Trees');
ylabel('Out-of-Bag Classification Error');
title('Out-of-Bag Error vs. Number of Grown Trees - PCA');

% Display the out-of-bag classification error of each individual tree
numTrees = rfModelPCA.NumTrees;
individualErrors = zeros(numTrees, 1);

for i = 1:numTrees
    % Get the out-of-bag predictions for the i-th tree
    [~, scores] = oobPredict(rfModelPCA, 'Trees', i);
    oobSamples = find(rfModelPCA.OOBIndices(:, i));  % Indices of OOB samples for the i-th tree
    oobPredictions = scores(oobSamples, :);       % Predictions for OOB samples
    [~, oobPredictedLabels] = max(oobPredictions, [], 2); % Convert scores to labels
    oobPredictedLabels = categorical(oobPredictedLabels, 1:length(categories(filteredLabels)), categories(filteredLabels));
    
    % Calculate OOB error for the i-th tree
    individualErrors(i) = sum(oobPredictedLabels ~= filteredLabels(oobSamples)) / length(oobSamples);
end

figure;
plot(1:numTrees, individualErrors, 'o-');
xlabel('Tree Index');
ylabel('Out-of-Bag Classification Error');
title('OOB Error of Individual Trees - PCA');

% Display predictor importance
predictorImportance = rfModelPCA.OOBPermutedPredictorDeltaError;
figure;
bar(predictorImportance);
xlabel('Predictor Index');
ylabel('Predictor Importance');
title('Predictor Importance Estimates - PCA');

% Display the structure of at least three trees
view(rfModelPCA.Trees{1}, 'Mode', 'graph');
view(rfModelPCA.Trees{2}, 'Mode', 'graph');
view(rfModelPCA.Trees{3}, 'Mode', 'graph');

% Save the Random Forest model
save(fullfile(mainDir, 'rfModel_PCA.mat'), 'rfModelPCA');


%%%%%%%%% REDUCTION BY PREDICTOR IMPORTANCE %%%%%%%%%%%%%%%
% Identify the most informative features based on their importance
[~, sortedIndices] = sort(predictorImportance, 'descend');
mostInformativeFeatures = sortedIndices(1:10); % Select top 10 features

% Retrain the classifier using only the most informative features
trainFeaturesReduced = trainFeatures(:, mostInformativeFeatures);
testFeaturesReduced = testFeatures(:, mostInformativeFeatures);

% Reduced number of trees for retraining
reducedNTrees = 1000;

% Create the TreeBagger object with reduced number of trees
tic;
rfModelReduced = TreeBagger(reducedNTrees, trainFeaturesReduced, trainLabels, ...
    'OOBPredictorImportance', 'on', ...
    'InBagFraction', inBagFraction);

% Predict on the test data with the reduced model
predictedLabelsReduced = predict(rfModelReduced, testFeaturesReduced);
classificationTimeFull = toc;
predictedLabelsReduced = categorical(predictedLabelsReduced);

% Compute accuracy of the reduced model
accuracyReduced = sum(predictedLabelsReduced == testLabels) / numel(testLabels);
disp(['Accuracy on test set with reduced model: ' num2str(accuracyReduced * 100) '%']);
disp(['Classification time: ' num2str(classificationTimeFull) ' seconds']);

% Create confusion matrix for the reduced model
confMatReduced = confusionmat(testLabels, predictedLabelsReduced);

% Display confusion matrix for the reduced model
figure;
heatmap(categories(filteredLabels), categories(filteredLabels), confMatReduced);
xlabel('Predicted Labels');
ylabel('True Labels');
title('Confusion Matrix (Reduced Model)');

% Display out-of-bag error for the reduced model
oobErrorsReduced = oobError(rfModelReduced);
figure;
plot(oobErrorsReduced);
xlabel('Number of Grown Trees');
ylabel('Out-of-Bag Classification Error');
title('Out-of-Bag Error vs. Number of Grown Trees (Reduced Model)');

% Display predictor importance for the reduced model
predictorImportanceReduced = rfModelReduced.OOBPermutedPredictorDeltaError;
figure;
bar(predictorImportanceReduced);
xlabel('Predictor Index');
ylabel('Predictor Importance');
title('Predictor Importance Estimates (Reduced Model)');

% Save the reduced Random Forest model
save(fullfile(mainDir, 'rfModelReduced.mat'), 'rfModelReduced');
