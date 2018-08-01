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

allqnn = {};
allPreds = {};
timing = [];
timingBase = [];

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

    out_est = localTempoEstimation(score{1}, expertDB, gettempo(performance));
    out_est(:,5) = deal(inputFolders(i));
    onsetEst = vertcat(out_est{:,1});
    
    allPreds = vertcat(allPreds, out_est);
        
    % midi generation for perceptual test
    predmid = vertcat(predqnn{:,1});
    predmid(predmid(:,5) < 1,5) = deal(1);
    predmid = predmid(:,1:7);
    predmid(:,1:2) = score{1}(:,1:2);
    predmid(:,6:7) = onsetEst(:,6:7);
    perfmid = performance(:,1:7);
    %perfmid(:,1:2) = performance(:,11:12);
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

    trivialSqErr = abs(performance(:,5) - score{2});
    nnSqErr = abs(nnDyn - performance(:,5));
    qnnSqErr = abs(qnnDyn - performance(:,5));
    knnSqErr = abs(knnDyn - performance(:,5));
    wknnSqErr = abs(wknnDyn - performance(:,5));

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

    timing = [timing; abs(predmid(:,7) - perfmid(:,7))];
    timingBase = [timingBase; abs(deadpmid(:,7) - perfmid(:,7))];
    
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

    qnn(:,7) = feats(:,1);
    allqnn = [allqnn; qnn];
end
dda = vertcat(da{1,:});
ddb = vertcat(db{1,:});
ddc = vertcat(dc{1,:});



for i = 1:1:size(allPreds,1)
    allPreds{i,6} = size(allPreds{i,1},1);
end

scoresSamples = ...
    {11, [1,2]; 12, [3,4,5,6]; 13,[7,8,9]; ...
     21, [10,11;5,5]; 22,[11,12;6,12]; 23,[11;6]; ...
     31, [13,14;11,15]; 32,[14,15,16,17,18;3,10,8,2,8]; 33,[18,19;2,14]; ...
     41, [20,21,22,23]; 42,[24,25,26,27,28,29]; 43,[30,31,32,33]; ...
     51, [35,36,37]; 52,[38]; 53,[39,40,41,42]; ...
     61, [43,44,45]; 62,[46,47,48,49;7,6,4,3]; 63,[49,50,51,52;3,3,5,6]; ...
     71, [53,54]; 72,[55,56]; 73,[57,58]; ...
     81, [59,60,61,62,63,64,65]; 82,[66,67,68;16,4,8]; 83, [68:73;1,2,5,2,4,16];
     };

for i = 1:size(scoresSamples,1)
    range = scoresSamples{i,2}(1,:);
    weights = [];
    if size(scoresSamples{i,2},1) > 1
        weights = scoresSamples{i,2}(2,:);
    else
        weights = [allPreds{scoresSamples{i,2}(1,:),6}];
    end
    scoresSamples{i,3} = [allPreds{range,2}]*weights'/sum(weights);
end



figure
boxplot([(base1), (nn1), (qnn1), (knn1)], ...
    'Notch', 'off', 'Labels', ...
    {'Baseline (mechanical)','kNN (exact, k=1)', 'kNN (parabola, k=1)','kNN (k=3)'});
title('Distribution of note level errors in dynamics predictions vs. baseline');
ylabel('Error in predicted note velocity (1-127)');

%figure
%boxplot([abs(base2), nn2, qnn2, knn2, wknn2], ...
%    'Notch', 'on', 'Labels', ...
%    {'Baseline','1-NN (exact)', '1-NN (parabola)','k-NN','weighted k-NN'});
%title('Distribution of note level errors in dynamics contour predictions vs. Baseline');
%ylabel('Error in predicted loudness z-score');

%figure
%boxplot([abs(base3), nn3, knn3, wknn3], ...
%    'Notch', 'on', 'Labels', ...
%    {'Baseline','1-NN','k-NN','weighted k-NN'});
%title('Distribution of phrase level errors in dynamic salience predictions vs. trivial approach');
%ylabel('Error in predicted loudness (????)');

for i = 1:size(allqnn,1)
    allqnn{i,8} = size(allqnn{i,1},1);
end
allqnn = allqnn([allqnn{:,8}] > 3,:);
lowdist = allqnn([allqnn{:,6}] < quantile([allqnn{:,6}], 0.3),:);
hidist = allqnn([allqnn{:,6}] > quantile([allqnn{:,6}], 0.7),:);
hidistDyn = vertcat(hidist{:,1});
hidistDyn = hidistDyn(:,5); 
hidistDyn(:,2:14) = vertcat(hidist{:,7});
hidistDyn = hidistDyn(:,[1 6]);
hidistDyn(:,3) = abs(hidistDyn(:,1) - hidistDyn(:,2));
hidistMAE = mean(hidistDyn(:,3));
lodistDyn = vertcat(lowdist{:,1});
lodistDyn = lodistDyn(:,5); 
lodistDyn(:,2:14) = vertcat(lowdist{:,7});
lodistDyn = lodistDyn(:,[1 6]);
lodistDyn(:,3) = abs(lodistDyn(:,1) - lodistDyn(:,2));
lodistMAE = mean(lodistDyn(:,3));