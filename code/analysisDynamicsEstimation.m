%% analysis of dynamics estimation

clear
%% Load the required libraries

addpath(genpath('miditoolbox'));

%% create expert database
trainingSet = {'beethoven4_4E1', 0}; %; 'bachManual', 0; 'meditacion', 0; 'borodin2_1', 0; 'haydn', 0};
expertDB = createExpertDB(trainingSet, 0);
s = size(expertDB,1);

if isunix(), sep = '/'; else sep = '\'; end

%%%%% TODO: compute overall lvl and dynamic range of test/xval sets under
%%%%% the premise that the model only has to guess the relative salience
%%%%% and contours.

% load score and calculate expressive features from performance
testSet = {'beethoven4_3', 0; 'beethoven4_1', 0};
performance = [];
score = {};
feats = {};
for i = 1:size(testSet, 1)
    f = createExpertDB(testSet(i,:), 0);
    p = vertcat(f{:,1});
    score{i} = p(:,1:7);
    performance = [performance; p]; %#ok<AGROW>
    feats = [feats; f]; %#ok<AGROW>
end

% load a cross-validation set
[XVscore, XVperf] = exprFeat(['..' sep 'data' sep 'bachManual' sep 'input'], 0, 0);

% estimate dynamics for cross-validation set and choose parameters
for i = 1:6
    for j = 1:5
        pred = dynamicsEstimation(XVscore, expertDB, 'w-knn', i+5, 10^(j-9));
        errnotes = vertcat(pred{:,1});
        errnotes = errnotes(:,5) - XVperf(:,5);
        sqerr(i,j) = errnotes'*errnotes;
    end
end

[aux, min_i] = min(sqerr);
[~, min_j] = min(min_i);
k = min_i(min_j)+5;
c = 10^(min_j-9);

%% estimate dynamics for target score

predwknn = {};
predknn = {};
predqnn = {};
prednn = {};
for i = 1:size(testSet,1)
    [wknn, knn, qnn, nn] = dynamicsEstimation(score{i}, expertDB, 'all', k, c);
    predwknn = [predwknn; wknn]; %#ok<AGROW>
    predknn = [predknn; knn]; %#ok<AGROW>
    predqnn = [predqnn; qnn]; %#ok<AGROW>
    prednn = [prednn; nn]; %#ok<AGROW>
end

clear folderList training wknn knn qnn nn aux min_i min_j i j pred errnotes f p

%% test quality of output according to original performance

% 1. Note-level analyses

% calculate overall mean level for comparison
trivialLevel = mean([expertDB{:,4}]);

% generate predicted dynamics curve from concatenation of all segments in
% predictions and calculate mean squared error
nnDyn = vertcat(prednn{:,1});
nnDyn = nnDyn(:,5);
qnnDyn = vertcat(predqnn{:,1});
qnnDyn = qnnDyn(:,5);
knnDyn = vertcat(predknn{:,1});
knnDyn = knnDyn(:,5);
wknnDyn = vertcat(predwknn{:,1});
wknnDyn = wknnDyn(:,5);

trivialSqErr = (performance(:,5) - trivialLevel).^2;
nnSqErr = (nnDyn - performance(:,5)).^2;
qnnSqErr = (qnnDyn - performance(:,5)).^2;
knnSqErr = (knnDyn - performance(:,5)).^2;
wknnSqErr = (wknnDyn - performance(:,5)).^2;

figure
boxplot([sqrt(trivialSqErr), sqrt(nnSqErr), sqrt(qnnSqErr), sqrt(knnSqErr), sqrt(wknnSqErr)], ...
    'Notch', 'on', 'Labels', ...
    {'Trivial (mean loudness)','1-NN (exact)', '1-NN (parabola)','k-NN','weighted k-NN'});
title('Distribution of note level errors in dynamics predictions vs. trivial approach');
ylabel('Error in predicted loudness (dBFS)');

performedContour = [];
for i = 1:size(feats,1)
    performedContour = vertcat(performedContour, feats{i,1}(:,5) - feats{i,4}); %#ok<AGROW>
end

nnContourErr = abs(vertcat(prednn{:,4}) - performedContour);
qnnContourErr = abs(vertcat(predqnn{:,4}) - performedContour);
knnContourErr = abs(vertcat(predknn{:,4}) - performedContour);
wknnContourErr = abs(vertcat(predwknn{:,4}) - performedContour);

figure
boxplot([abs(performedContour), nnContourErr, qnnContourErr, knnContourErr, wknnContourErr], ...
    'Notch', 'on', 'Labels', ...
    {'Baseline','1-NN (exact)', '1-NN (parabola)','k-NN','weighted k-NN'});
title('Distribution of note level errors in dynamics contour predictions vs. Baseline');
ylabel('Error in predicted loudness (dBFS)');

clear performedContour nnContourErr qnnContourErr wknnContourErr

performedSalience = vertcat(feats{:,6});
nnErr = abs(vertcat(prednn{:,5}) - performedSalience);
knnErr = abs(vertcat(predknn{:,5}) - performedSalience);
wknnErr = abs(vertcat(predwknn{:,5}) - performedSalience);

figure
boxplot([abs(performedSalience), nnErr, knnErr, wknnErr], ...
    'Notch', 'on', 'Labels', ...
    {'Trivial (mean loudness)','1-NN','k-NN','weighted k-NN'});
title('Distribution of phrase level errors in dynamic salience predictions vs. trivial approach');
ylabel('Error in predicted loudness (dBFS)');

clear performedSalience nnErr knnErr wknnErr
