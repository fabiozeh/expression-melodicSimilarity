%% File generation for perceptual test

clear
%% Load the required libraries

addpath(genpath('miditoolbox'));

%% List input folders desired
inputFolders = {'ev01'; 'ev02'; 'ev03'; 'ev04'; 'ev05'; 'ev06'; 'ev07'; 'ev09'};

%% Statistics tables
base1 = [];
nn1 = [];
qnn1 = [];
knn1 = [];
wknn1 = [];
base2 = [];
nn2 = [];
qnn2 = [];
knn2 = [];
wknn2 = [];
base3 = [];
nn3 = [];
qnn3 = [];
knn3 = [];
wknn3 = [];

%% Generate models for all folders with a leave-one-out training set
for i = 1:length(inputFolders)
    % create expert database
    trainingSet = inputFolders;
    trainingSet(i,:) = []; % leave this out
    trainingSet(:,2) = {0};
    expertDB = createExpertDB(trainingSet, 0, 1);
    s = size(expertDB,1);

    if isunix(), sep = '/'; else sep = '\'; end

    % load score and calculate expressive features from performance
    testSet = {inputFolders{i}, 0};
    score = cell(1, 3);
    [feats, fData] = createExpertDB(testSet, 0, 0);
    performance = vertcat(feats{:,1});
    score{1} = performance(:,1:7);
    score{1}(:,5) = 80; % neutralize velocities for fairness (TODO test if makes difference)
    % compute piece overall mean dynamics
    score{2} = fData(1,1);
    % compute piece dynamic range
    score{3} = fData(1,2);

    % estimate dynamics for target score
    k = 3;
    c = 1e-6;
    predwknn = {};
    predknn = {};
    predqnn = {};
    prednn = {};
    
    [wknn, knn, qnn, nn] = dynamicsEstimation(score{1}, 64, 100, expertDB, 'all', k, c);
    predwknn = [predwknn; wknn]; %#ok<AGROW>
    predknn = [predknn; knn]; %#ok<AGROW>
    predqnn = [predqnn; qnn]; %#ok<AGROW>
    prednn = [prednn; nn]; %#ok<AGROW>

    score{1}(:,6:7) = score{1}(:,1:2) ./ mean(performance(:,13)) * 60; % tempo adjustment
    out_est = onsetDevEstimation(score{1}, fData(1,3), fData(1,4), expertDB);
    onsetEst = vertcat(out_est{:,1});
    
    % midi generation for perceptual test
    predmid = vertcat(predqnn{:,1});
    predmid = predmid(:,1:7);
    predmid(:,1:2) = score{1}(:,1:2);
    predmid(:,6:7) = onsetEst(:,6:7);
    perfmid = performance(:,1:7);
    perfmid(:,1:2) = performance(:,11:12);
    perfmid(:,6:7) = performance(:,9:10);
    
    deadpmid = score{1};
    deadpmid(:,5) = score{2};
    
    da(i) = num2cell(perfmid, [1 2]);
    db(i) = num2cell(predmid, [1 2]);
    dc(i) = num2cell(deadpmid, [1 2]);
    writemidi(deadpmid, [inputFolders{i} '_deadpan.mid']);
    writemidi(perfmid, [inputFolders{i} '_performed.mid']);
    writemidi(predmid, [inputFolders{i} '_predicted.mid']);

    % data analysis
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

    trivialSqErr = (performance(:,5) - score{2}).^2;
    nnSqErr = (nnDyn - performance(:,5)).^2;
    qnnSqErr = (qnnDyn - performance(:,5)).^2;
    knnSqErr = (knnDyn - performance(:,5)).^2;
    wknnSqErr = (wknnDyn - performance(:,5)).^2;

    performedContour = [];
    for i = 1:size(feats,1)
        performedContour = vertcat(performedContour, feats{i,10}); %#ok<AGROW>
    end

    nnContourErr = abs(vertcat(prednn{:,5}) - performedContour);
    qnnContourErr = abs(vertcat(predqnn{:,5}) - performedContour);
    knnContourErr = abs(vertcat(predknn{:,5}) - performedContour);
    wknnContourErr = abs(vertcat(predwknn{:,5}) - performedContour);

    performedSalience = vertcat(feats{:,8});
    nnSalErr = abs(vertcat(prednn{:,3}) - performedSalience);
    qnnSalErr = abs(vertcat(predqnn{:,3}) - performedSalience);
    knnSalErr = abs(vertcat(predknn{:,3}) - performedSalience);
    wknnSalErr = abs(vertcat(predwknn{:,3}) - performedSalience);

%     piano = quantile(performance(:,5),0.33);
%     forte = quantile(performance(:,5),0.66);
% 
%     perfp = performance(:,5) < piano;
%     %perfn = performance(:,5) >= quantile(performance(:,5),0.33) && ...
%     %    performance(:,5) <= quantile(performance(:,5),0.66);
%     perff = performance(:,5) > forte;
% 
%     perfQ = -1.*perfp + perff;
% 
%     knnDQ = -1.*(knnDyn < piano) + (knnDyn > forte);
% 
%     pctRight = sum(perfQ == knnDQ)./size(perfQ,1);
    
    
    base1 = [base1; trivialSqErr]; %#ok<AGROW>
    nn1 = [nn1; nnSqErr];%#ok<AGROW>
    qnn1 = [qnn1; qnnSqErr];%#ok<AGROW>
    knn1 = [knn1; knnSqErr];%#ok<AGROW>
    wknn1 = [wknn1; wknnSqErr];%#ok<AGROW>

    base2 = [base2; performedContour];%#ok<AGROW>
    nn2 = [nn2; nnContourErr];%#ok<AGROW>
    qnn2 = [qnn2; qnnContourErr];%#ok<AGROW>
    knn2 = [knn2; knnContourErr];%#ok<AGROW>
    wknn2 = [wknn2; wknnContourErr];%#ok<AGROW>

    base3 = [base3; performedSalience];%#ok<AGROW>
    nn3 = [nn3; nnSalErr];%#ok<AGROW>
    qnn3 = [qnn3; qnnSalErr];%#ok<AGROW>
    knn3 = [knn3; knnSalErr];%#ok<AGROW>
    wknn3 = [wknn3; wknnSalErr];%#ok<AGROW>

end
dda = vertcat(da{1,:});
ddb = vertcat(db{1,:});
ddc = vertcat(dc{1,:});

figure
boxplot([sqrt(base1), sqrt(nn1), sqrt(qnn1), sqrt(knn1), sqrt(wknn1)], ...
    'Notch', 'on', 'Labels', ...
    {'Baseline (mechanical)','1-NN (exact)', '1-NN (parabola)','k-NN','weighted k-NN'});
title('Distribution of note level errors in dynamics predictions vs. baseline');
ylabel('Error in predicted loudness (1-127)');

figure
boxplot([abs(base2), nn2, qnn2, knn2, wknn2], ...
    'Notch', 'on', 'Labels', ...
    {'Baseline','1-NN (exact)', '1-NN (parabola)','k-NN','weighted k-NN'});
title('Distribution of note level errors in dynamics contour predictions vs. Baseline');
ylabel('Error in predicted loudness z-score');

figure
boxplot([abs(base3), nn3, knn3, wknn3], ...
    'Notch', 'on', 'Labels', ...
    {'Baseline','1-NN','k-NN','weighted k-NN'});
title('Distribution of phrase level errors in dynamic salience predictions vs. trivial approach');
ylabel('Error in predicted loudness (????)');