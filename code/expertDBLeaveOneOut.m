%% Perform the prediction and analysis
    
%% Load the required libraries

addpath(genpath('miditoolbox'));

%% Load Expert Database
createExpertDB

%% calculate expertDB leave-one-out analysis

% calculate the matrix of segment similarities
scores = Inf(size(expertDB, 1));
for jj = 1:size(scores,2)
    for ii = 1:size(scores,1)
        if ii < jj
            scores(ii, jj) = scores(jj, ii); % it is simmetric
        elseif ii > jj % for equal indices, score = 0, so skip.
            [H, tbk] = dtwSig2(expertDB{jj,1}, expertDB{ii,1}, 0, 1, 0, 1, 'no');
            scores(ii, jj) = H(tbk(end,1),tbk(end,2)); %sum(h(sub2ind(size(h),tbk(:,1)+1,tbk(:,2)+1)))
        end
    end
end

clear ii jj H tbk

%% Compute dynamic estimations for each segment

for ii = 1:size(scores,1)
    [expertDB{ii,4}, ind] = min(scores(ii,:)); % for the lowest value of each row (could be column)
    perfStartInd = expertDB{ind,2};
    % take mean level in segment
    mL = mean(expertDB{ind,1}(:,5));
    expertDB{ii,5} = scaleInterpolate(size(expertDB{ii,1},1), expertDB{ind,1}(:,5) - mL);
    % expertDB{:,6} --> offset from previous section
    if perfStartInd > 1
        expertDB{ii,6} = mL - mean(expertDB{ind-1,1}(:,5));
    else
        expertDB{ii,6} = 0;
    end
    expertDB{ii,7} = mL; % segment mean level
    expertDB{ii,8} = mean(cellfun(@(x) mean(x(:,5)), expertDB{expertDB{:,3}==expertDB{ind,3},1})); % piece mean level
end
% should interpolation be done in time domain instead of note domain?

% TODO
% more stable algorithm:
% if segment level is > 1 sd from mean, copy neighbor gradient
% else copy deviation from mean

% adjust mean level of segments
% the strategy is as follows:
% - two parameters are saved for each reference segment: the mean level of
% the piece it belongs to (meanVel = expertDB{:,8}) and the offset of the segment mean relative to
% the mean level of previous segment (os = expertDB{:,6}).
% - after determining all predictions, the final piece mean level is adjusted to
% be the mean(meanVel) taking all segments into account and the offset of each
% segment is corrected to os = os - mean(os).
meanOS = mean([segments{:,4}]);
offset = meanVel + segments{1,4} - meanOS;
segments{1,1}(:,5) = segments{1,1}(:,5) + offset;
for ii = 2:size(segments,1)
    offset = segments{ii,4} + offset - meanOS;
    segments{ii,1}(:,5) = segments{ii,1}(:,5) + offset;
end

clear ii ind perfStartInd perfEndInd meanOS mL meanVel offset

%% test quality of output according to original performance

% cc = [segmentCorrelation segmentMeanSqError segmentDistance]
cc(:,1) = arrayfun(@(x) corr(segments{x,1}(:,5), ...
    alignedperf(segments{x,2}:(segments{x,2}+size(segments{x,1},1)-1),5)), ...
    1:size(segments,1));
cc(:,2) = arrayfun(@(x) mean((segments{x,1}(:,5) - ...
    alignedperf(segments{x,2}:(segments{x,2}+size(segments{x,1},1)-1),5)).^2), ...
    1:size(segments,1));
cc(:,3) = arrayfun(@(x) segments{x,3}, 1:size(segments,1));

%ccEvolution = [maxDistance meanCorrelation meanSquareError numSegments]
% "for all segments with melodic distance d <= ..."
ccEvolution(:,1) = [0; 0.25; 0.5; 1; 2; 4; 8; 16; 32; 64];
% "the mean correlation between prediction and performance"
ccEvolution(:,2) = arrayfun(@(x)mean(cc(cc(:,3)<=x,1)), ccEvolution(:,1));
% "the mean squared error between prediction and performance"
ccEvolution(:,3) = arrayfun(@(x)mean(cc(cc(:,3)<=x,2)), ccEvolution(:,1));
% "the number of segments in this category (distance <= 1st col value)"
ccEvolution(:,4) = arrayfun(@(x)sum(cc(:,3)<=x), ccEvolution(:,1));


% calculate trivial (average-based) model for comparison
trivialModel = segments;
for ii = 1:size(trivialModel,1)
    perfStartInd = trivialModel{ii,2};
    perfEndInd = trivialModel{ii,2} + size(trivialModel{ii,1},1) - 1;
    v = mean([alignedperf(1:perfStartInd-1,5); alignedperf(perfEndInd+1:end,5)]);
    trivialModel{ii,1}(:,5) = v;
end
% trivial model predicts for every note in a segment, the mean velocity
% among notes in all other segments in the piece.
outputTrivial = cat(1, trivialModel{:,1});

trivMeanSqErr(1,:) = arrayfun(@(x) mean((trivialModel{x,1}(:,5) - ...
    alignedperf(trivialModel{x,2}:(trivialModel{x,2}+size(trivialModel{x,1},1)-1),5)).^2), ...
    1:size(trivialModel,1));

% ccT_Ev(i) = the mean squared error between trivial model and performance for all
% segments with distance d <= ccEvolution(i,1).
ccT_Ev(1,:) = arrayfun(@(x)mean(trivMeanSqErr(cc(:,3)<=x)), ccEvolution(:,1));

clear ii v perfStartInd perfEndInd trivialModel trivMeanSqErr


% This graph shows how mean squared error between prediction and
% performance within segments changes as the melodic distance of the match
% for each segment increases, whereas the same doesn't happen for the 
% trivial model which acts as a baseline. For distances below 2 we see that
% our model outperforms the trivial model in this criteria.
figure
hold on
plot(ccEvolution(:,1),ccEvolution(:,3));
plot(ccEvolution(:,1),ccT_Ev);
hold off

% This graph show the same fact as before, but instead of using mean
% squared error, we use mean correlation between predicted dynamics and
% performed dynamics within each segment. (obs: Here we can't compare to
% the trivial model because it generates the same value for all notes in
% the same segment and therefore no correlation value can be computed.)
%plot(ccEvolution(:,1),ccEvolution(:,2));

% the melodic distance of each note, for plotting against output curve
melDistance = Inf(size(outputmidi,1),1);
for ii = 1:size(segments,1)
    melDistance(segments{ii,2}:(segments{ii,2}+size(segments{ii,1},1)-1),1) = ...
        segments{ii,3} + zeros(size(segments{ii,1},1),1);
end
clear ii

% smooth functions for better visualization
truth = filter(0.1*ones(1,10),1,alignedperf(:,5));
prediction = filter(0.1*ones(1,10),1,outputmidi(:,5));
trivial = filter(0.1*ones(1,10),1,outputTrivial(:,5));
truth = truth(10:end); %discard the first elements because of filtering logic
prediction = prediction(10:end);
trivial = trivial(10:end);

errSmooth = (truth - prediction).^2;
errTrivial = (truth - trivial).^2;

% correlation between computed distance of a segment and the correlation
% value between prediction and performance of the same segment
corr_corrXmelDist = corr(cc(:,1),cc(:,3));
%plot(ccEvolution(:,2))


% correlation between the mean squared error between prediction and
% performance of a segment and the computed distance between score segment
% and reference segment.
corr_errXmelDist = corr(cc(:,2),cc(:,3));
%plot(ccEvolution(:,3))

%hold on
%plot(outputmidi(:,5));
%plot(alignedperf(:,5));
%hold off

figure
hold on
plot(prediction);
plot(truth);
hold off

figure
hold on
plot(errSmooth);
plot(errTrivial);
plot(melDistance(10:end));
hold off

%% reported data
x = cc(cc(:,3)<=2,1);
SEM = std(x)/sqrt(length(x));               % Standard Error
ts = tinv([0.025  0.975],length(x)-1);      % T-Score
PLMIN = ts*SEM;
CI = mean(x) + PLMIN;                       % Confidence Interval

x = cc(:,1);
SEM = std(x)/sqrt(length(x));               % Standard Error
ts = tinv([0.025  0.975],length(x)-1);      % T-Score
PLMINall = ts*SEM;
CIall = mean(x) + PLMINall;                 % Confidence Interval

%addpath('util');
%writeall(segments, [dataFolder sep 'analysis']);

clear x SEM ts PLMIN PLMINall dataFolder sep