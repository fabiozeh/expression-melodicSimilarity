%% generate all arffs with velocity, onset deviations and duration deviations

clear

%% Get desired input folder
if isunix(), sep = '/'; else sep = '\'; end
inputFolder = uigetdir(['.' sep], 'Select input data folder:');

addpath(genpath(['..' sep 'util' sep 'score2arff']));
addpath(genpath(['..' sep]));

[db, info] = expressionDataset(inputFolder);

if isunix()
    system('mkdir arff');
else
    system('md arff');
end

for i = 1:length(db)
    array2arff(db{i}, ['arff' sep int2str(i) 'test.arff'], info(:,1), info(:,2), info(:,3));
    dbi = db;
    dbi(i) = []; % delete this instance
    dbi = vertcat(dbi{:,:});
    array2arff(dbi, ['arff' sep int2str(i) 'train.arff'], info(:,1), info(:,2), info(:,3));
end

array2arff(vertcat(db{:,:}), ['arff' sep 'all.arff'], info(:,1), info(:,2), info(:,3));

