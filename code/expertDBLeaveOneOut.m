%% Perform the prediction and analysis
    
%% Load the required libraries

addpath(genpath('miditoolbox'));

%% Load Expert Database
createExpertDB
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
            scores(ii, jj) = H(tbk(end,1),tbk(end,2)); %sum(h(sub2ind(size(h),tbk(:,1)+1,tbk(:,2)+1)))
        end
    end
end

clear ii jj H tbk

% meanVel = { name of piece, mean piece velocity}
meanVel = expertDB([expertDB{:,2}] == 1, 3);
for ii = 1:size(meanVel,1)
    meanVel{ii,2} = mean([expertDB{strcmp(meanVel{ii,1}, expertDB(:,3)), 4}]);
end

% allocation
predictions{s,9} = [];

for ii = 1:size(scores,1)
    [predictions{ii,2}, ind] = min(scores(ii,:)); % for the lowest value of each row (could be column)
    predictions{ii,1} = ind;
    % base level prediction as linear interpolation of matching segment levels
    predictions{ii,3} = scaleInterpolate(size(expertDB{ii,1},1), expertDB{ind,1}(:,5) - expertDB{ind,4});
    % matching segment offset from previous
    predictions{ii,4} = expertDB{ind,5};
    % matching segment mean level
    predictions{ii,5} = expertDB{ind,4};
    % matching piece mean level
    predictions{ii,6} = meanVel{strcmp(expertDB{ind,3}, meanVel(:,1)), 2};
    % normalized melodic distance
    predictions{ii,7} = predictions{ii,2}/size(expertDB{ii,1},1);
end
% should interpolation be done in time domain instead of note domain?

% adjust mean level of segments

% the algorith below simply makes the segment's mean level deviate from the
% piece mean level by the same amount as the matching (reference) segment
% deviates from its original piece mean level.
%outMean = mean([predictions{:,6}]);
%for ii = 1:size(predictions,1)
%    predictions{ii,8} = predictions{ii,3} + outMean + predictions{ii,5} - predictions{ii,6};
%end

% more stable algorithm:
% if ref. segment level is > 1 sd from mean, copy neighbor gradient
% else copy deviation from mean
outMean = mean([predictions{:,6}]);
for ii = 1:size(predictions,1)
    if expertDB{predictions{ii,1},6} > 1
        predictions{ii,9} = outMean + predictions{ii,4};
    else
        predictions{ii,9} = outMean + predictions{ii,5} - predictions{ii,6};
    end
    predictions{ii,8} = predictions{ii,3} + predictions{ii,9};
end

clear ind meanVel

%% test quality of output according to original performance

% 1. Note-level analyses

% calculate overall mean level for comparison
trivialLevel = mean([expertDB{:,4}]);

% generate performed dynamics curve from concatenation of all segments in DB
performedDyn = vertcat(expertDB{:,1});
performedDyn = performedDyn(:,5);

% generate predicted dynamics curve from concatenation of all segments in DB
predictedDyn = vertcat(predictions{:,8});

trivialSqErr = (performedDyn - trivialLevel).^2;
predictedSqErr = (predictedDyn - performedDyn).^2;

figure
boxplot([sqrt(predictedSqErr),sqrt(trivialSqErr)],'Notch', 'on', 'Labels',{'Model','Trivial (average value)'});
title('Distribution of note level errors in dynamics predictions vs. trivial approach');
ylabel('Error in predicted loudness (0-127)');

clear trivialSqErr predictedSqErr

performedInSegmentDyn = [];
for ii = 1:s
    performedInSegmentDyn = vertcat(performedInSegmentDyn, expertDB{ii,1}(:,5) - expertDB{ii,4}); %#ok<AGROW>
end
predictedInSegmentDyn = vertcat(predictions{:,3});

figure
boxplot([abs(predictedInSegmentDyn - performedInSegmentDyn),abs(performedInSegmentDyn)],'Notch', 'on', 'Labels',{'Model','Baseline'});
title('Distribution of note level errors in short-range dynamics predictions vs. Baseline');
ylabel('Error in predicted loudness (0-127)');

clear predictedInSegmentDyn performedInSegmentDyn

predictedInSegmentDynMSE = zeros(s, 1);
segmentRMS = zeros(s,1);
for ii = 1:s
    predictedInSegmentDynMSE(ii) = mean((expertDB{ii,1}(:,5) - expertDB{ii,4} - predictions{ii,3}).^2);
    segmentRMS(ii) = rms(expertDB{ii,1}(:,5) - expertDB{ii,4});
end

figure
boxplot([sqrt(predictedInSegmentDynMSE),segmentRMS], ...
    'Notch', 'on', 'Labels',{'Model','Baseline'});
title('Distribution of mean note level errors in short-range dynamics predictions vs. Baseline');
ylabel('Error in predicted loudness (0-127)');

% correlation between intra segment level predictions MSE and melodic distance
corr_inSegDynMseXmelDist = corr(predictedInSegmentDynMSE, vertcat(predictions{:,6}));

clear predictedInSegmentDynMSE segmentRMS segmentDiffRMS

% 2. Segment-level analyses

% 2.1 best-match prediction
performedDyn = vertcat(expertDB{:,4});
predictedDyn = vertcat(predictions{:,9});

figure
boxplot([abs(predictedDyn - performedDyn), abs(performedDyn - trivialLevel)], ...
    'Notch', 'on', 'Labels',{'Model','Trivial (average-based)'});
title('Distribution of phrase level errors in dynamics predictions vs. trivial approach');
ylabel('Error in predicted loudness (0-127)');

% correlation between segment mean level prediction errors and melodic distance
corr_errSegLevelXmelDist = corr(abs(predictedDyn - performedDyn), vertcat(predictions{:,7}));

% 2.2 weighted sum prediction

% weights are proportional to the reciprocal of the melodic distance
w = 1./(scores+1e-6);
% some normalization
w = w./(ones(size(w))*diag(std(w)));
w = w.^2;

wsPred = [expertDB{:,4}]*w./(sum(w));

figure
boxplot([abs(wsPred' - performedDyn), abs(performedDyn - trivialLevel)], ...
    'Notch', 'on', 'Labels',{'Model','Trivial (average-based)'});
title('Distribution of phrase level errors in weighted-sum dynamics predictions vs. trivial approach');
ylabel('Error in predicted loudness (0-127)');

% 3. Evidence of melodic distance validity as a metric

% 3.1 For segment mean dynamics level
sq_err_lvl = [];
for ii = 1:s
    l = size(expertDB{ii,1}, 1);
    for jj = ii+1:s
        sq_err_lvl = [sq_err_lvl; ii, jj, scores(ii,jj)/l, abs(expertDB{ii, 6} - expertDB{jj, 6})]; %#ok<AGROW>
    end
end
corr_segZXallDist = corr(sq_err_lvl(:,3),sq_err_lvl(:,4));

for ii = 1:s
    x = w(:,ii);
    y = abs(expertDB{ii, 6} - [expertDB{:, 6}])';
    x(ii) = [];
    y(ii) = [];
    corrOneSegErrXweight(ii) = corr(x, y);
end

% 3.2 For note-level dynamics

% weighted-sum prediction based on quadratic regression
for ii = 1:s
    quadCoef(ii,:) = polyfit(scaleInterpolate(size(expertDB{ii,1},1), [0;10]), ...
        expertDB{ii,1}(:,5) - expertDB{ii,4}, 2);
end

wCoef = w * quadCoef ./ (sum(w,2)*ones(1,3));
for ii = 1:s
    segx = scaleInterpolate(size(expertDB{ii,1},1), [0;10]);
    predictions{ii,10} = (segx.^2).*wCoef(ii,1) + segx.*wCoef(ii,2) + wCoef(ii,3) + wsPred(ii);
end

for ii = 1:s
    segx = scaleInterpolate(size(expertDB{ii,1},1), [0;10]);
    c = quadCoef(predictions{ii,1},:);
    predictions{ii,11} = segx.^2.*c(1) + segx.*c(2) + c(3) + outMean;
end
clear segx c

polyfitMetrics = zeros(s,4);
for ii = 1:s
    polyfitMetrics(ii,1) = corr(expertDB{ii,1}(:,5), predictions{ii,10});
    polyfitMetrics(ii,2) = mean((expertDB{ii,1}(:,5) - predictions{ii,10}).^2);
    polyfitMetrics(ii,3) = corr(expertDB{ii,1}(:,5), predictions{ii,11});
    polyfitMetrics(ii,4) = mean((expertDB{ii,1}(:,5) - predictions{ii,11}).^2);
end

% Missing: sort segments by melodic distance and calculate moving average
% of errors to improve upon metrics evolution below. 

% metrics = [segmentCorrelation segmentMeanSqError normalizedSegmentDistance segmentNumberOfNotes corrP-val]
metrics = zeros(s,4);
for ii = 1:s
    metrics(ii,1) = corr(expertDB{ii,1}(:,5), predictions{ii,8});
    metrics(ii,2) = mean((expertDB{ii,1}(:,5) - predictions{ii,8}).^2);
    metrics(ii,3) = predictions{ii,7}; % change to 2 for non-normalized
    metrics(ii,4) = size(expertDB{ii,1},1);
end

% metricsEvolution = [maxDistance meanCorrelation meanSquareError numSegments]
metricsEvolution(:,1) = [-0.05; 0; 0.1; 0.2; 0.4; 0.8; 1.6; 3.2; 6.4; 1000];
metricsEvolution(:,2:4) = deal(NaN); % preallocation
for ii = 2:size(metricsEvolution,1)
    % "for all segments with melodic distance ... < d <= ..."
    f = metrics(:,3) <= metricsEvolution(ii,1) & metrics(:,3) > metricsEvolution(ii-1,1);
    % "the mean correlation between prediction and performance"
    metricsEvolution(ii,2) = mean(metrics(f,1));
    % "the mean squared error between prediction and performance"
    metricsEvolution(ii,3) = mean(metrics(f,2));
    % "the number of segments in this category (distance <= 1st col value)"
    metricsEvolution(ii,4) = sum(f);
end
metricsEvolution = metricsEvolution(2:end,:);

metricsCumulative(:,1) = [0; 0.1; 0.2; 0.4; 0.8; 1.6; 3.2; 6.4; 1000]; % TODO change to quantiles
metricsCumulative(:,2:4) = deal(0);
for ii = 1:size(metricsCumulative,1)
    % "for all segments with melodic distance d <= ..."
    f = metrics(:,3) <= metricsCumulative(ii,1);
    % "the mean correlation between prediction and performance"
    metricsCumulative(ii,2) = mean(metrics(f,1));
    % "the mean squared error between prediction and performance"
    metricsCumulative(ii,3) = mean(metrics(f,2));
    % "the number of segments in this category (distance <= 1st col value)"
    metricsCumulative(ii,4) = sum(f);
end

metricsCumuAll(:,1) = [0; 0.1; 0.2; 0.4; 0.8; 1.6; 3.2; 6.4; 12.8; 50; 100; 200; 1000];
metricsCumuAll(:,2:4) = deal(0);
for ii = 1:size(metricsCumuAll,1)
    f = sq_err_lvl(:,3) <= metricsCumuAll(ii,1);
    if sum(f) ~= 0
        metricsCumuAll(ii,2) = mean(sq_err_lvl(f,4));
        metricsCumuAll(ii,3) = corr(sq_err_lvl(f,3),sq_err_lvl(f,4));
        metricsCumuAll(ii,4) = sum(f);
    end
end
    
notesEvolution(:,1) = [0; 2; 4; 6; 8; 10; 12; 16; 100];
notesEvolution(:,2:4) = deal(NaN); % preallocation
for ii = 2:size(notesEvolution,1)
    % "for all segments with number of notes ... < n <= ..."
    f = metrics(:,4) <= notesEvolution(ii,1) & metrics(:,4) > notesEvolution(ii-1,1);
    % "the mean correlation between prediction and performance"
    notesEvolution(ii,2) = mean(metrics(f,1));
    % "the mean squared error between prediction and performance"
    notesEvolution(ii,3) = mean(metrics(f,2));
    % "the number of segments in this category"
    notesEvolution(ii,4) = sum(f);
end
notesEvolution = notesEvolution(2:end,:);

trivialSegmentMse = arrayfun(@(x) mean((expertDB{x,1}(:,5) - trivialLevel).^2), 1:s);

% mseTrivDistCumulative(i) = the mean squared error between trivial model and performance for all
% segments with distance d <= ccEvolution(i,1).
mseTrivDistCumulative(1,:) = arrayfun(@(x)mean(trivialSegmentMse(metrics(:,3)<=x)), metricsEvolution(:,1));

clear f ii trivialLevel trivialSegmentMse

% This graph shows how mean squared error between prediction and
% performance within segments changes as the melodic distance of the match
% for each segment increases.
figure
hold on
plot(metricsEvolution(:,3));
plot(mseTrivDistCumulative);
hold off

% This graph shows how mean correlation between prediction vs. performance
% loudness varies as we take into account segments for which the matching
% sample is increasingly distant melodically (thus worse).
figure
plot(metricsCumulative(:,2)); % why are all correlations so high? plot some and try everyone against everyone

% This calculates the note-level correlations among all segment pairs
corrNoteLvlAll = zeros(s);
for ii = 1:s
    for jj = 1:s
        corrNoteLvlAll(ii,jj) = ...
            corr(expertDB{ii,1}(:,5) - expertDB{ii,4}, ...
            scaleInterpolate(size(expertDB{ii,1},1), expertDB{jj,1}(:,5) - expertDB{jj,4}));
    end
end
% It can be seen that the best possible correlation for each segment is
% most times far superior to the one chosen by 1-NN with melodic distance
% metric (3 in 56 coincide. Expected coincidences in random bernoulli trials
% should be 56/55). Still, mean correlation of chosen segments is much
% higher than mean correlation overall.

c_cnlXweight = zeros(s,1);
for ii = 1:s
    x = w(:,ii);
    y = corrNoteLvlAll(:,ii);
    x(ii) = [];
    y(ii) = [];
    c_cnlXweight(ii) = corr(x, y);
end
[areNoteLevelCorrAndWeightCorrelated, ~, c_cnlXweight_ci, ~] = ttest(c_cnlXweight);

clear x y ii

% the melodic distance of each note, for plotting against output curve
melDistance = Inf(size(performedDyn,1),1);
for ii = 1:s
    melDistance(expertDB{ii,2}:(expertDB{ii,2}+size(expertDB{ii,1},1)-1),1) = ...
        deal(expertDB{ii,6});
end
clear ii

% correlation between computed distance of a segment and the correlation
% value between prediction and performance of the same segment.
corr_corrXmelDist = corr(metrics(:,1),metrics(:,3));

% correlation between the mean squared error between prediction and
% performance of a segment and the computed distance between score segment
% and reference segment.
corr_errXmelDist = corr(metrics(:,2),metrics(:,3));

% correlation between prediction vs. performance correlation and number of
% notes in segment. (desired = 0)
corr_corrXnotes = corr(metrics(:,1),metrics(:,4));

% correlation between prediction vs. performance mean squared error and
% number of notes in segment. 
corr_errXnotes = corr(metrics(:,2),metrics(:,4));

x = metrics(metrics(:,3)<=0.4,1);
SEM = std(x)/sqrt(length(x));               % Standard Error
ts = tinv([0.025  0.975],length(x)-1);      % T-Score
PLMIN = ts*SEM;
CI = mean(x) + PLMIN;                       % Confidence Interval

x = metrics(:,1);
SEM = std(x)/sqrt(length(x));               % Standard Error
ts = tinv([0.025  0.975],length(x)-1);      % T-Score
PLMINall = ts*SEM;
CIall = mean(x) + PLMINall;                 % Confidence Interval

%addpath('util');
%writeall(segments, [dataFolder sep 'analysis']);

clear x SEM ts PLMIN PLMINall dataFolder sep s