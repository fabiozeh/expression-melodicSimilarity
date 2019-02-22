%% plots the learning rate by progressively increasing the dataset.

clear
if isunix(), sep = '/'; else, sep = '\'; end

%% Load the required libraries
addpath(genpath(['..' sep 'miditoolbox']));
addpath(genpath(['..' sep]));

%% Get desired input folder
inputFolder = uigetdir(['.' sep], 'Select input data folder:');
%inputFolder = '../../data/feb19larger'; 

%% Load training set and compute expressive features
fullDB = createExpertDB(inputFolder, 0, 0);
S = size(fullDB,1);

%% randomize samples
rng(997); % setting seed
perm = randperm(S);
fullDB = fullDB(perm,:);
[~, sortvec] = sort(perm); % for putting it back in order

%% set number of cross-validation folds and learning curve datapoints
folds = 10;
datapoints = 10;
exponentialIncrease = 1;
firstDatapointSize = 5*folds;

%% choose values of k for kNN
k = [5 6 7 8];

%% Statistics tables
mae = cell(datapoints,length(k),5); % 3 values of k, 5 measurements
mae_train = cell(datapoints,length(k),5);
pred_all = cell(datapoints,2);

ratio = exp(log(S/(firstDatapointSize))/datapoints);
for run = 1:datapoints

    if exponentialIncrease
        s = firstDatapointSize*ratio.^run;
    else
        s = max(folds, floor(S*run/datapoints)); %#ok
    end
    thisDB = fullDB(1:s,:);
    foldsize = max(1,floor(s/folds));
    largerfolds = s - folds*foldsize;

    i = 1;
    predgroup = [];
    traingroup = [];
    foldsize = foldsize + 1;
    while i < s

        if i == largerfolds*foldsize + 1
            foldsize = foldsize - 1;
        end

        xval = thisDB(i:min(s,i+foldsize-1),:);
        train = vertcat(thisDB(1:i-1,:), thisDB(min(s,i+foldsize):end,:));

        i = i + foldsize;

        % estimate dynamics for cross-validation samples. 
        % Mean vel (arg 2) = mean of xval pieces mean dynamics
        knn = dynamicsEstimation(xval, mean([xval{:,12}]), 30, train, 'knn', k);
        
        % estimate dynamics for training samples
        trainset = train(randperm(size(xval,1)),:);
        train_knn = dynamicsEstimation(trainset, mean([trainset{:,12}]), 30, train, 'knn', k);

        % data analysis
        
        % generate predicted dynamics curve from concatenation of all segments
        dynPreds = cell(length(k),2);
        for kind = 1:(length(k))
            velvals = vertcat(knn{:,1,kind});
            dynPreds{kind,1} = velvals(:,5);
        
            velvals = vertcat(train_knn{:,1,kind});
            dynPreds{kind,2} = velvals(:,5);
        end
        
        xval_groundtruth = vertcat(xval{:,1});
        xval_groundtruth = dbfs2vel_sqrt(xval_groundtruth(:,5)); % set velocities in midi vals

        train_groundtruth = vertcat(trainset{:,1});
        train_groundtruth = dbfs2vel_sqrt(train_groundtruth(:,5));
        
        % compute mean absolute errors
        for kind = 1:length(k)
            mae{run,kind,1} = [mae{run,kind,1}; abs(dynPreds{kind,1} - xval_groundtruth)]; %output velocity values
            mae{run,kind,2} = [mae{run,kind,2}; abs(vertcat(xval{:,8}) - vertcat(knn{:,3,kind}))]; %alpha
            mae{run,kind,3} = [mae{run,kind,3}; abs(vertcat(xval{:,9}) - vertcat(knn{:,4,kind}))]; %beta
            mae{run,kind,4} = [mae{run,kind,4}; abs(vertcat(xval{:,16}) - vertcat(knn{:,6,kind}))]; %gamma coefs
            mae{run,kind,5} = [mae{run,kind,5}; vertcat(knn{:,7,kind})]; % mean distance
            mae_train{run,kind,1} = [mae_train{run,kind,1}; abs(dynPreds{kind,2} - train_groundtruth)];
            mae_train{run,kind,2} = [mae_train{run,kind,2}; abs(vertcat(trainset{:,8}) - vertcat(train_knn{:,3,kind}))]; %alpha
            mae_train{run,kind,3} = [mae_train{run,kind,3}; abs(vertcat(trainset{:,9}) - vertcat(train_knn{:,4,kind}))]; %beta
            mae_train{run,kind,4} = [mae_train{run,kind,4}; abs(vertcat(trainset{:,16}) - vertcat(train_knn{:,6,kind}))]; %gamma coefs
            mae_train{run,kind,5} = [mae_train{run,kind,5}; vertcat(train_knn{:,7,kind})]; % mean distance
        end
        predgroup = [predgroup; knn]; %#ok<AGROW>
        traingroup = [traingroup; train_knn]; %#ok<AGROW>
    end
    pred_all{run,1} = predgroup;
    pred_all{run,2} = traingroup;
end

means = zeros(datapoints,3,length(k));
fstQ = zeros(datapoints,3,length(k));
trdQ = zeros(datapoints,3,length(k));
for i = 1:datapoints
    for kind = 1:length(k)
        means(i,1,kind) = size(mae{i,kind,2},1);
        means(i,2,kind) = median(mae{i,kind,1});
        means(i,3,kind) = median(mae_train{i,kind,1});
        fstQ(i,1,kind) = quantile(mae{i,kind,1},0.25);
        fstQ(i,2,kind) = quantile(mae_train{i,kind,1},0.25);
        trdQ(i,1,kind) = quantile(mae{i,kind,1},0.75);
        trdQ(i,2,kind) = quantile(mae_train{i,kind,1},0.75);
    end
end

plot_kind = 2;

baseline = zeros(size(means,1),1);
for i = 1:size(means,1)
    s = means(i,1,plot_kind);
    fulldyn = vertcat(fullDB{1:s,1});
    fulldyn = dbfs2vel_sqrt(fulldyn(:,5));
    meandyn = 0;
    s = floor(size(fulldyn,1)/10);
    for j = 1:10
        meandyn = meandyn + mean(fulldyn((j-1)*s+1:j*s,1))./10;
    end
    baseline(i) = median(abs(meandyn - fulldyn));
end

%plotting only k = 3
figure1 = figure;
hold on

xfill = [means(:,1,plot_kind); flipud(means(:,1,plot_kind))];
yfill1 = [fstQ(:,1,plot_kind); flipud(trdQ(:,1,plot_kind))];
fill(xfill, yfill1, [0.25 0.25 0.92], 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
yfill2 = [fstQ(:,2,plot_kind); flipud(trdQ(:,2,plot_kind))];
fill(xfill, yfill2, [0.92 0.25 0.25], 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(means(:,1,plot_kind),fstQ(:,1,plot_kind),means(:,1,plot_kind),fstQ(:,2,plot_kind), 'LineStyle',':', 'Color', [0.25 0.25 0.92], 'HandleVisibility', 'off');
plot(means(:,1,plot_kind),trdQ(:,1,plot_kind),means(:,1,plot_kind),trdQ(:,2,plot_kind), 'LineStyle',':', 'Color', [0.92 0.25 0.25], 'HandleVisibility', 'off');

line([68 68],[0.2 20], 'Color', [0 0 0], 'HandleVisibility', 'off'); %DS1
line([192 192],[2 20], 'Color', [0 0 0], 'HandleVisibility', 'off'); %DS2
line([2706 2706],[0.2 20], 'Color', [0 0 0], 'HandleVisibility', 'off'); %DS3

plot(means(:,1,plot_kind), means(:,2,plot_kind), ...
    'LineWidth', 2, 'Color', [0.25 0.25 0.92], 'DisplayName','cross-validation set (median)');

plot(means(:,1,plot_kind), means(:,3,plot_kind), ...
    'LineWidth', 2, 'Color', [0.92 0.25 0.25], 'DisplayName','training set (median)');

plot(means(:,1,plot_kind), baseline(:,1), ...
    'LineWidth', 2, 'Color', [0.25 0.9 0.9], 'DisplayName','deadpan (median)');

title('Evolution of MAE with dataset size');
ylabel('Error in predicted note velocity (1-127)');
xlabel('Number of instances (motifs)');

axes1 = gca;
xlim(axes1,[0 2800]);
box(axes1,'on');
legend(axes1,'show');

% Create textbox
annotation(figure1,'textbox',...
    [0.436323383084577 0.731277533039651 0.0300945273631841 0.0418502202643172],...
    'Color',[1 0 0],...
    'String','Q3',...
    'LineWidth',1,...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.433835820895522 0.233480176211454 0.0300945273631841 0.0418502202643172],...
    'Color',[1 0 0],...
    'String','Q1',...
    'LineWidth',1,...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.389059701492537 0.861233480176213 0.0300945273631841 0.0418502202643172],...
    'Color',[0 0 1],...
    'String','Q3',...
    'LineWidth',1,...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.38781592039801 0.312775330396476 0.0300945273631841 0.0418502202643172],...
    'Color',[0 0 1],...
    'String','Q1',...
    'LineWidth',1,...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.153985074626866 0.12114537444934 0.0761144278606966 0.0418502202643171],...
    'String','DS1 size',...
    'LineWidth',1,...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.187567164179104 0.200440528634362 0.0761144278606966 0.0418502202643171],...
    'String','DS2 size',...
    'LineWidth',1,...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.806970149253731 0.151982378854626 0.0761144278606966 0.0418502202643171],...
    'String','DS3 size',...
    'LineWidth',1,...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor','none');
