%% Perform the prediction and analysis

clear
%% Load the required libraries

addpath(genpath('miditoolbox'));

%% Load Expert Database
dataSet = {'beethoven4_4E1', 0}; %; 'beethoven4_3', 0; 'beethoven4_1', 0; 'bachManual', 0}%; 'meditacion', 0; 'borodin2_1', 0; 'haydn', 0};
[expertDB, dbData] = createExpertDB(dataSet, 0);

% delete the segment with a single note
expertDB(end,:) = []; 


s = size(expertDB,1);

% exclude small segments (2 notes)
% expertDB = expertDB(arrayfun(@(x) size(expertDB{x,1},1) > 2, 1:s),:);

%% calculate expertDB leave-one-out analysis

% calculate the matrix of segment similarities
scores = Inf(s);
for jj = 1:s
    for ii = 1:s
        if ii < jj
            scores(ii, jj) = scores(jj, ii); % it is simmetric
        elseif ii > jj % for equal indices, leave it out, so skip.
            [H, tbk] = dtwSig2(expertDB{jj,1}, expertDB{ii,1}, 0, 1, 0, 1, 'no');
            scores(ii, jj) = H(tbk(end,1),tbk(end,2));
        end
    end
end

% compute the quadratic coefficients that approximate the contour curve for
% each expertDB segment

quadCoef = zeros(s,3);
expertDBxs = cell(s,1);
for i = 1:s
    x = expertDB{i,1};
    x0 = x(1,1);
    x1 = x(end,1)+x(end,2);
    x = x(:,1)+0.5.*x(:,2);
    expertDBxs{i} = lin_interpolation(x, 0, x0, 10, x1);
    quadCoef(i,:) = polyfit(expertDBxs{i}, ...
        expertDB{i,10},2);
end

clear ii jj H tbk

% allocation
predictions{s,13} = [];
for ii = 1:s
    [predictions{ii,2}, ind] = min(scores(ii,:)); % for the lowest value of each row (could be column)
    predictions{ii,1} = ind;
    % nearest-neighbor alpha
    predictions{ii,3} = expertDB{ind,8};
    % nearest-neighbor beta
    predictions{ii,4} = expertDB{ind,9};
    % nearest-neighbor gamma
    predictions{ii,5} = lin_approx(expertDBxs{ii}, expertDBxs{ind}, expertDB{ind,10});
    % nearest-neighbor gamma by quadratic approximation
    predictions{ii,6} = quadCoef(ind,1).*expertDBxs{ii}.^2 + quadCoef(ind,2).* ...
        expertDBxs{ii} + quadCoef(ind,3);
    % 1-nn prediction
    predictions{ii,7} = dbData(1) + dbData(2).*(predictions{ii,3} + ...
        predictions{ii,4}.*predictions{ii,5});
    % 1-nn prediction with quadratic approximation
    predictions{ii,8} = dbData(1) + dbData(2).*(predictions{ii,3} + ...
        predictions{ii,4}.*predictions{ii,6});
    % normalized melodic distance
    predictions{ii,9} = predictions{ii,2}/size(expertDB{ii,1},1);
end

clear ind

% weights are proportional to the reciprocal of the melodic distance
w = 1./(scores+1e-9);
% some normalization
w = w./(diag(std(w,0,2))*ones(size(w)));
%w = w.^2;
% kNN
k = 3;
w2 = zeros(s);
for ii = 1:k
    [v, ind] = max(w,[],2);
    w(sub2ind([s,s],1:s,ind')) = 0;
    w2(sub2ind([s,s],1:s,ind')) = 1; %(1 for equal weights, v for w-kNN)
end
w = w2;
clear k ind w2

wCoef = w * quadCoef ./ (sum(w,2)*ones(1,3));

% k-NN Alpha
predictions(:,10) = num2cell(w*vertcat(expertDB{:,8})./sum(w,2));
% k-NN Beta
predictions(:,11) = num2cell(w*vertcat(expertDB{:,9})./sum(w,2));

for ii = 1:s
    % k-NN gamma
    predictions{ii,12} = (expertDBxs{ii}.^2).*wCoef(ii,1) + ...
        expertDBxs{ii}.*wCoef(ii,2) + wCoef(ii,3);
    % k-NN prediction
    predictions{ii,13} = dbData(1) + dbData(2).*(predictions{ii,10} + ...
        predictions{ii,11}.*predictions{ii,12});
end

%% test quality of output according to original performance

% 1. Note-level analyses

% calculate overall mean level for comparison
trivialLevel = dbData(1);

% generate performed dynamics curve from concatenation of all segments in DB
performedDyn = vertcat(expertDB{:,1});
performedDyn = performedDyn(:,5);

% generate predicted dynamics curve from concatenation of all segments in
% predictions and calculate mean squared error

trivialErr = sqrt((performedDyn - trivialLevel).^2);
nnErr = sqrt((vertcat(predictions{:,7}) - performedDyn).^2);
nnqErr = sqrt((vertcat(predictions{:,8}) - performedDyn).^2);
knnErr = sqrt((vertcat(predictions{:,13}) - performedDyn).^2);

figure
boxplot([trivialErr, nnErr, nnqErr, knnErr], ...
    'Notch', 'on', 'Labels', ...
    {'Baseline (mechanical)','1-NN', '1-NN (parabola)','k-NN'});
title('Distribution of note level errors in dynamics predictions');
ylabel('Error in predicted dynamics (0-127)');

% p-value reported
[~, knnPval, ~, ~] = ttest(trivialErr, knnErr);

performedGamma = vertcat(expertDB{:,10});

nnGammaErr = abs(vertcat(predictions{:,5}) - performedGamma);
nnqGammaErr = abs(vertcat(predictions{:,6}) - performedGamma);
knnGammaErr = abs(vertcat(predictions{:,12}) - performedGamma);

figure
boxplot([abs(performedGamma), nnGammaErr, nnqGammaErr, knnGammaErr], ...
    'Notch', 'on', 'Labels', ...
    {'Baseline','1-NN', '1-NN (parabola)','k-NN'});
title('Distribution of note level errors in dynamics contour predictions vs. Baseline');
ylabel('Error in predicted loudness (????)');

% 2. Segment-level analyses

performedAlpha = vertcat(expertDB{:,8});
nnAlphaErr = abs(vertcat(predictions{:,3}) - performedAlpha);
knnAlphaErr = abs(vertcat(predictions{:,10}) - performedAlpha);

figure
boxplot([abs(performedAlpha), nnAlphaErr, knnAlphaErr], ...
    'Notch', 'on', 'Labels', ...
    {'Baseline','1-NN', 'k-NN'});
title('Distribution of segment level errors in alpha predictions vs. Baseline');
ylabel('Error in predicted loudness (????)');

tableResults = [ median(trivialErr), median(nnErr), median(nnqErr), median(knnErr), knnPval ];

% separating phrases with closest nn from those with farthest nn.
maeseg = Inf(s,1);
corrseg = Inf(s,1);
for i=1:s
    perfseg = expertDB{i,1}(:,5);
    predseg = predictions{i,7};
    maeseg(i) = mean(abs(predseg - perfseg));
    corrseg(i) = corr(perfseg, predseg);
end
[nnScoreOrd, sortorder] = sort(vertcat(predictions{:,2}));
nnMaeSegOrd = maeseg(sortorder);

% correlation between melodic distance and 1-nn mean absolute error
%corr(nnScoreOrd,nnMaeSegOrd)

% correlation between melodic distance and 1-nn predicted segment 
% correlation w/ performance
%corr(nnScoreOrd,nnMaeSegOrd)

half = floor(s/2);
if half < s/2, last = s - 1; else last = s; end

% mean absolute error in 1-nn for closest phrases in range 0 - 1.
mean(nnMaeSegOrd(1:half)) / 127 

% mean absolute error in 1-nn for farthest phrases in range 0 - 1.
mean(nnMaeSegOrd(half+1:last)) / 127

% 3. Evidence of melodic distance validity as a metric

% 3.1 For segment mean dynamics level
% sq_err_lvl = [];
% for ii = 1:s
%     l = size(expertDB{ii,1}, 1);
%     for jj = ii+1:s
%         sq_err_lvl = [sq_err_lvl; ii, jj, scores(ii,jj)/l, abs(expertDB{ii, 6} - expertDB{jj, 6})]; %#ok<AGROW>
%     end
% end
% corr_segZXallDist = corr(sq_err_lvl(:,3),sq_err_lvl(:,4));
% 
% for ii = 1:s
%     x = w(:,ii);
%     y = abs(expertDB{ii, 6} - [expertDB{:, 6}])';
%     x(ii) = [];
%     y(ii) = [];
%     corrOneSegErrXweight(ii) = corr(x, y);
% end
% 
% % 3.2 For note-level dynamics
% 
% polyfitMetrics = zeros(s,4);
% for ii = 1:s
%     polyfitMetrics(ii,1) = corr(expertDB{ii,1}(:,5), predictions{ii,13});
%     polyfitMetrics(ii,2) = mean((expertDB{ii,1}(:,5) - predictions{ii,13}).^2);
%     polyfitMetrics(ii,3) = corr(expertDB{ii,1}(:,5), predictions{ii,14});
%     polyfitMetrics(ii,4) = mean((expertDB{ii,1}(:,5) - predictions{ii,14}).^2);
% end
% 
% % Missing: sort segments by melodic distance and calculate moving average
% % of errors to improve upon metrics evolution below. 
% 
% % metrics = [segmentCorrelation segmentMeanSqError normalizedSegmentDistance segmentNumberOfNotes corrP-val]
% metrics = zeros(s,4);
% for ii = 1:s
%     metrics(ii,1) = corr(expertDB{ii,1}(:,5), predictions{ii,8});
%     metrics(ii,2) = mean((expertDB{ii,1}(:,5) - predictions{ii,8}).^2);
%     metrics(ii,3) = predictions{ii,7}; % change to 2 for non-normalized
%     metrics(ii,4) = size(expertDB{ii,1},1);
% end
% 
% % metricsEvolution = [maxDistance meanCorrelation meanSquareError numSegments]
% metricsEvolution(:,1) = [-0.05; 0; 0.1; 0.2; 0.4; 0.8; 1.6; 3.2; 6.4; 1000];
% metricsEvolution(:,2:4) = deal(NaN); % preallocation
% for ii = 2:size(metricsEvolution,1)
%     % "for all segments with melodic distance ... < d <= ..."
%     f = metrics(:,3) <= metricsEvolution(ii,1) & metrics(:,3) > metricsEvolution(ii-1,1);
%     % "the mean correlation between prediction and performance"
%     metricsEvolution(ii,2) = mean(metrics(f,1));
%     % "the mean squared error between prediction and performance"
%     metricsEvolution(ii,3) = mean(metrics(f,2));
%     % "the number of segments in this category (distance <= 1st col value)"
%     metricsEvolution(ii,4) = sum(f);
% end
% metricsEvolution = metricsEvolution(2:end,:);
% 
% metricsCumulative(:,1) = [0; 0.1; 0.2; 0.4; 0.8; 1.6; 3.2; 6.4; 1000]; % TODO change to quantiles
% metricsCumulative(:,2:4) = deal(0);
% for ii = 1:size(metricsCumulative,1)
%     % "for all segments with melodic distance d <= ..."
%     f = metrics(:,3) <= metricsCumulative(ii,1);
%     % "the mean correlation between prediction and performance"
%     metricsCumulative(ii,2) = mean(metrics(f,1));
%     % "the mean squared error between prediction and performance"
%     metricsCumulative(ii,3) = mean(metrics(f,2));
%     % "the number of segments in this category (distance <= 1st col value)"
%     metricsCumulative(ii,4) = sum(f);
% end
% 
% metricsCumuAll(:,1) = [0; 0.1; 0.2; 0.4; 0.8; 1.6; 3.2; 6.4; 12.8; 50; 100; 200; 1000];
% metricsCumuAll(:,2:4) = deal(0);
% for ii = 1:size(metricsCumuAll,1)
%     f = sq_err_lvl(:,3) <= metricsCumuAll(ii,1);
%     if sum(f) ~= 0
%         metricsCumuAll(ii,2) = mean(sq_err_lvl(f,4));
%         metricsCumuAll(ii,3) = corr(sq_err_lvl(f,3),sq_err_lvl(f,4));
%         metricsCumuAll(ii,4) = sum(f);
%     end
% end
%     
% notesEvolution(:,1) = [0; 2; 4; 6; 8; 10; 12; 16; 100];
% notesEvolution(:,2:4) = deal(NaN); % preallocation
% for ii = 2:size(notesEvolution,1)
%     % "for all segments with number of notes ... < n <= ..."
%     f = metrics(:,4) <= notesEvolution(ii,1) & metrics(:,4) > notesEvolution(ii-1,1);
%     % "the mean correlation between prediction and performance"
%     notesEvolution(ii,2) = mean(metrics(f,1));
%     % "the mean squared error between prediction and performance"
%     notesEvolution(ii,3) = mean(metrics(f,2));
%     % "the number of segments in this category"
%     notesEvolution(ii,4) = sum(f);
% end
% notesEvolution = notesEvolution(2:end,:);
% 
% trivialSegmentMse = arrayfun(@(x) mean((expertDB{x,1}(:,5) - trivialLevel).^2), 1:s);
% 
% % mseTrivDistCumulative(i) = the mean squared error between trivial model and performance for all
% % segments with distance d <= ccEvolution(i,1).
% mseTrivDistCumulative(1,:) = arrayfun(@(x)mean(trivialSegmentMse(metrics(:,3)<=x)), metricsEvolution(:,1));
% 
% clear f ii trivialLevel trivialSegmentMse
% 
% % This graph shows how mean correlation between prediction vs. performance
% % loudness varies as we take into account segments for which the matching
% % sample is increasingly distant melodically (thus worse).
% figure
% plot(metricsCumulative(:,2)); % why are all correlations so high? plot some and try everyone against everyone
% 
% % This calculates the note-level correlations among all segment pairs
% corrNoteLvlAll = zeros(s);
% for ii = 1:s
%     for jj = 1:s
%         corrNoteLvlAll(ii,jj) = ...
%             corr(expertDB{ii,1}(:,5) - expertDB{ii,4}, ...
%             scaleInterpolate(size(expertDB{ii,1},1), expertDB{jj,1}(:,5) - expertDB{jj,4}));
%     end
% end
% % It can be seen that the best possible correlation for each segment is
% % most times far superior to the one chosen by 1-NN with melodic distance
% % metric (3 in 56 coincide. Expected coincidences in random bernoulli trials
% % should be 56/55). Still, mean correlation of chosen segments is much
% % higher than mean correlation overall.
% 
% c_cnlXweight = zeros(s,1);
% low_score = 0.4;
% y1 = [];
% for ii = 1:s
%     x = scores(ii,:);
%     y = corrNoteLvlAll(ii,:);
%     y1 = [y1 corrNoteLvlAll(ii, scores(ii,:) <= low_score)];
%     x(ii) = [];
%     y(ii) = [];
%     c_cnlXweight(ii) = corr(x', y');
% end
% [areNoteLevelCorrAndWeightCorrelated, ~, c_cnlXweight_ci, ~] = ttest(c_cnlXweight);
% [~, ~, nlc_highW_ci, ~] = ttest(y1');
% 
% clear x y ii jj
% 
% % the melodic distance of each note, for plotting against output curve
% melDistance = Inf(size(performedDyn,1),1);
% for ii = 1:s
%     melDistance(expertDB{ii,2}:(expertDB{ii,2}+size(expertDB{ii,1},1)-1),1) = ...
%         deal(expertDB{ii,6});
% end
% clear ii
% 
% % correlation between computed distance of a segment and the correlation
% % value between prediction and performance of the same segment.
% corr_corrXmelDist = corr(metrics(:,1),metrics(:,3));
% 
% % correlation between the mean squared error between prediction and
% % performance of a segment and the computed distance between score segment
% % and reference segment.
% corr_errXmelDist = corr(metrics(:,2),metrics(:,3));
% 
% % correlation between prediction vs. performance correlation and number of
% % notes in segment. (desired = 0)
% corr_corrXnotes = corr(metrics(:,1),metrics(:,4));
% 
% % correlation between prediction vs. performance mean squared error and
% % number of notes in segment. 
% corr_errXnotes = corr(metrics(:,2),metrics(:,4));
% 
% x = metrics(metrics(:,3)<=0.4,1);
% SEM = std(x)/sqrt(length(x));               % Standard Error
% ts = tinv([0.025  0.975],length(x)-1);      % T-Score
% PLMIN = ts*SEM;
% CI = mean(x) + PLMIN;                       % Confidence Interval
% 
% x = metrics(:,1);
% SEM = std(x)/sqrt(length(x));               % Standard Error
% ts = tinv([0.025  0.975],length(x)-1);      % T-Score
% PLMINall = ts*SEM;
% CIall = mean(x) + PLMINall;                 % Confidence Interval

%addpath('util');
%writeall(segments, [dataFolder sep 'analysis']);



clear v x0 x1 x last i ii half SEM ts PLMIN PLMINall dataFolder sep s