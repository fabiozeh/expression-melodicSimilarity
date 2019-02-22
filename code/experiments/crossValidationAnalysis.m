%% Phrase-level modeling cross-validation

clear
if isunix(), sep = '/'; else, sep = '\'; end

%% Load the required libraries
addpath(genpath(['..' sep 'miditoolbox']));
addpath(genpath(['..' sep]));

%% Get desired input folder
inputFolder = uigetdir(['.' sep], 'Select input data folder:');

%% Statistics tables
base1 = [];
nn1 = [];
qnn1 = [];
knn1 = [];

%timing = [];
%timingBase = [];

%% Load training set and compute expressive features
expertDB = createExpertDB(inputFolder, 0, 0);
s = size(expertDB,1);

%% randomize samples
rng(997); % setting seed
perm = randperm(s);
expertDB = expertDB(perm,:);
[~, sortvec] = sort(perm); % for putting it back in order

folds = 10;
foldsize = max(1,floor(s/folds));
largerfolds = s - folds*foldsize;

% outputs
da = {};
db = {};
dc = {};

i = 1;
foldsize = foldsize + 1;
while i < s
    
    if i == largerfolds*foldsize + 1
        foldsize = foldsize - 1;
    end
    
    xval = expertDB(i:min(s,i+foldsize-1),:);
    train = vertcat(expertDB(1:i-1,:), expertDB(min(s,i+foldsize):end,:));
    
    i = i + foldsize;

    performance = vertcat(xval{:,1});
    performance = dbfs2vel_sqrt(performance); % set velocities in midi vals
    score{1} = performance(:,1:7);
    % mean of xval pieces mean dynamics
    score{2} = mean([xval{:,12}]);

    % estimate dynamics for target score
    k = 2:7;
    c = 1e-6;
    
    [~, knn, qnn, nn] = dynamicsEstimation(xval, score{2}, 30, train, 'all', k, c);
    
    %out_est = localTempoEstimation(score, train, gettempo(performance));
    %onsetEst = vertcat(out_est{:,1});
    
    % midi generation
    predmid = vertcat(qnn{:,1});
    predmid = predmid(:,1:7);
    predmid(:,1:2) = score{1}(:,1:2);
    %predmid(:,6:7) = onsetEst(:,6:7);
    perfmid = performance(:,1:7);
    perfmid(:,6:7) = performance(:,9:10);
    
    deadpmid = score{1};
    deadpmid(:,5) = dbfs2vel_sqrt(score{2});
    
    da = vertcat(da, num2cell(perfmid, [1 2])); %#ok<AGROW>
    db = vertcat(db, num2cell(predmid, [1 2])); %#ok<AGROW>
    dc = vertcat(dc, num2cell(deadpmid, [1 2])); %#ok<AGROW>
    
    % data analysis
    % generate predicted dynamics curve from concatenation of all segments in
    % predictions and calculate mean squared error
    nnDyn = vertcat(nn{:,1});
    nnDyn = nnDyn(:,5);
    qnnDyn = vertcat(qnn{:,1});
    qnnDyn = qnnDyn(:,5);
    knnDyn = zeros(size(qnnDyn,1),length(k));
    for kind = 1:length(k)
        knnDynAux = vertcat(knn{:,1,kind});
        knnDyn(:,kind) = knnDynAux(:,5);
    end
    clear knnDynAux
    
    trivialSqErr = abs(dbfs2vel_sqrt(score{2}) - performance(:,5));
    nnSqErr = abs(nnDyn - performance(:,5));
    qnnSqErr = abs(qnnDyn - performance(:,5));
    knnSqErr = abs(knnDyn - performance(:,5));

%    timing = [timing; abs(predmid(:,7) - perfmid(:,7))];
%    timingBase = [timingBase; abs(deadpmid(:,7) - perfmid(:,7))];
    
    base1 = [base1; trivialSqErr]; %#ok<AGROW>
    nn1 = [nn1; nnSqErr];%#ok<AGROW>
    qnn1 = [qnn1; qnnSqErr];%#ok<AGROW>
    knn1 = [knn1; knnSqErr];%#ok<AGROW>
end

figure
boxplot([(base1), (nn1), (qnn1), (knn1(:,2))], ...
    'Notch', 'off', 'Labels', ...
    {'Baseline (deadpan)','kNN (exact, k=1)', 'kNN (parabola, k=1)','kNN (k=3)'});
title('Distribution of errors in dynamics predictions');
ylabel('Error in predicted note velocity (1-127)');
