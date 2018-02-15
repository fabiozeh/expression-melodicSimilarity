%% File generation for perceptual test

clear
%% Load the required libraries

addpath(genpath('miditoolbox'));

%% List input folders desired
%inputFolders = {'haydn'; 'beethoven4_1'; 'beethoven4_4E1'; 'borodin2_1'; 
    %'beethoven4_3'; 'meditacion'; 'bachManual'; 'ev01'; 'ev09'};

    inputFolders = {'ev01'; 'ev02'; 'ev03'; 'ev06'; 'ev09'};
    
%% Generate models for all folders with a leave-one-out training set
for i = 1:length(inputFolders)
    % create expert database
    trainingSet = inputFolders;
    trainingSet(i,:) = []; % leave this out
    trainingSet(:,2) = {0};
    expertDB = createExpertDB(trainingSet, 0);
    s = size(expertDB,1);

    if isunix(), sep = '/'; else sep = '\'; end

    % load score and calculate expressive features from performance
    testSet = {inputFolders{i}, 0};
    performance = [];
    score = cell(size(testSet, 1), 3);
    feats = {};
    for j = 1:size(testSet, 1)
        [f, fData] = createExpertDB(testSet(j,:), 0);
        p = vertcat(f{:,1});
        score{j,1} = p(:,1:7);
        % compute piece overall mean dynamics
        score{j,2} = fData(1,1);
        % compute piece dynamic range
        score{j,3} = fData(1,2);
        performance = [performance; p]; %#ok<AGROW>
        feats = [feats; f]; %#ok<AGROW>
    end

    % estimate dynamics for target score
    k = 3;
    c = 1e-6;
    predknn = {};
    for j = 1:size(testSet,1)
        knn = dynamicsEstimation(score{j,1}, score{j,2}, score{j,3}, expertDB, 'qnn', k, c);
        predknn = [predknn; knn]; %#ok<AGROW>
    end

    predmid = vertcat(predknn{:,1});
    predmid = predmid(:,1:7);
    perfmid = performance(:,1:7);
    deadpmid = perfmid;
    deadpmid(:,5) = 80;
    
    writemidi(deadpmid, [inputFolders{i} '_deadpan.mid']);
    writemidi(perfmid, [inputFolders{i} '_performed.mid']);
    writemidi(predmid, [inputFolders{i} '_predicted.mid']);
end
