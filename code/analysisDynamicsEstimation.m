%% analysis of dynamics estimation

clear
%% Load the required libraries

addpath(genpath('miditoolbox'));

%% create expert database
trainingSet = {'beethoven4_4E1', 0}; %'meditacion', 0; 'borodin2_1', 0; 'haydn', 0}; %'bachManual', 0; 
expertDB = createExpertDB(trainingSet, 0);
s = size(expertDB,1);

if isunix(), sep = '/'; else sep = '\'; end

%%%%% TODO: compute overall lvl and dynamic range of test/xval sets under
%%%%% the premise that the model only has to guess the relative salience
%%%%% and contours.

% load score and calculate expressive features from performance
testSet = {'beethoven4_3', 0}%; 'beethoven4_1', 0};
performance = [];
score = cell(size(testSet, 1), 3);
feats = {};
for i = 1:size(testSet, 1)
    [f, fData] = createExpertDB(testSet(i,:), 0);
    p = vertcat(f{:,1});
    score{i,1} = p(:,1:7);
    % compute piece overall mean dynamics
    score{i,2} = fData(1,1);
    % compute piece dynamic range
    score{i,3} = fData(1,2);
    performance = [performance; p]; %#ok<AGROW>
    feats = [feats; f]; %#ok<AGROW>
end

% load a cross-validation set
% [XVscore, XVperf] = exprFeat(['..' sep 'data' sep 'bachManual' sep 'input'], 0, 0);
% 
% % estimate dynamics for cross-validation set and choose parameters
% cvals = [1e-6 1e-3 1 100 1e4];
% kvals = [3 5 7 9 11];
% for i = 1:5
%     for j = 1:5
%         pred = dynamicsEstimation(XVscore, expertDB, 'w-knn', kvals(i), cvals(j));
%         errnotes = vertcat(pred{:,1});
%         errnotes = (errnotes(:,5) - XVperf(:,5)).^2;
%         sqerr(i,j) = median(errnotes);
%     end
% end
% 
% [aux, min_i] = min(sqerr);
% [~, min_j] = min(aux);
% k = kvals(min_i(min_j));
% c = cvals(min_j);

k = 5;
c = 1e-6;

%% estimate dynamics for target score

predwknn = {};
predknn = {};
predqnn = {};
prednn = {};
for i = 1:size(testSet,1)
    [wknn, knn, qnn, nn] = dynamicsEstimation(score{i,1}, score{i,2}, score{i,3}, expertDB, 'all', k, c);
    predwknn = [predwknn; wknn]; %#ok<AGROW>
    predknn = [predknn; knn]; %#ok<AGROW>
    predqnn = [predqnn; qnn]; %#ok<AGROW>
    prednn = [prednn; nn]; %#ok<AGROW>
end

clear folderList training wknn knn qnn nn min_i min_j aux cvals kvals i j pred errnotes f fData p

%% test quality of output according to original performance

% 1. Note-level analyses

% calculate overall mean level for comparison
trivialLevel = mean(performance(:,5));

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
    {'Baseline (mechanical)','1-NN (exact)', '1-NN (parabola)','k-NN','weighted k-NN'});
title('Distribution of note level errors in dynamics predictions vs. baseline');
ylabel('Error in predicted loudness (0-127)');

performedContour = [];
for i = 1:size(feats,1)
    performedContour = vertcat(performedContour, feats{i,10}); %#ok<AGROW>
end

nnContourErr = abs(vertcat(prednn{:,5}) - performedContour);
qnnContourErr = abs(vertcat(predqnn{:,5}) - performedContour);
knnContourErr = abs(vertcat(predknn{:,5}) - performedContour);
wknnContourErr = abs(vertcat(predwknn{:,5}) - performedContour);

figure
boxplot([abs(performedContour), nnContourErr, qnnContourErr, knnContourErr, wknnContourErr], ...
    'Notch', 'on', 'Labels', ...
    {'Baseline','1-NN (exact)', '1-NN (parabola)','k-NN','weighted k-NN'});
title('Distribution of note level errors in dynamics contour predictions vs. Baseline');
ylabel('Error in predicted loudness z-score');

clear performedContour nnContourErr qnnContourErr wknnContourErr

performedSalience = vertcat(feats{:,8});
nnErr = abs(vertcat(prednn{:,3}) - performedSalience);
knnErr = abs(vertcat(predknn{:,3}) - performedSalience);
wknnErr = abs(vertcat(predwknn{:,3}) - performedSalience);

figure
boxplot([abs(performedSalience), nnErr, knnErr, wknnErr], ...
    'Notch', 'on', 'Labels', ...
    {'Baseline','1-NN','k-NN','weighted k-NN'});
title('Distribution of phrase level errors in dynamic salience predictions vs. trivial approach');
ylabel('Error in predicted loudness (????)');

clear performedSalience nnErr knnErr wknnErr

piano = quantile(performance(:,5),0.33);
forte = quantile(performance(:,5),0.66);

perfp = performance(:,5) < piano;
%perfn = performance(:,5) >= quantile(performance(:,5),0.33) && ...
%    performance(:,5) <= quantile(performance(:,5),0.66);
perff = performance(:,5) > forte;

perfQ = -1.*perfp + perff;

knnDQ = -1.*(knnDyn < piano) + (knnDyn > forte);

pctRight = sum(perfQ == knnDQ)./size(perfQ,1);