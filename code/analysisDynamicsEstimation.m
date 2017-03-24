%% analysis of dynamics estimation

%% create expert database
createExpertDB

if isunix(), sep = '/'; else sep = '\'; end

folder = ['..' sep 'data' sep 'meditacion' sep 'input'];

%% Load the required libraries and files for score and performance

addpath(genpath('miditoolbox'));

% load score and calculate expressive features from performance
[score, performance] = exprFeat(folder, 0, 1);

%% estimate dynamics for target score

segments = dynamicsEstimation(score, expertDB);

% reconstruct predicted performance from segments
outputmidi = cat(1, segments{:,1});

%% test quality of output according to original performance

% cc = [segmentCorrelation segmentMeanSqError segmentDistance]
cc(:,1) = arrayfun(@(x) corr(segments{x,1}(:,5), ...
    performance(segments{x,2}:(segments{x,2}+size(segments{x,1},1)-1),5)), ...
    1:size(segments,1));
cc(:,2) = arrayfun(@(x) mean((segments{x,1}(:,5) - ...
    performance(segments{x,2}:(segments{x,2}+size(segments{x,1},1)-1),5)).^2), ...
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

mse = (outputmidi(:,5) - performance(:,5)).^2;
allsegs = vertcat(expertDB{:,1});
trivial = mean(allsegs(:,5));
trivialMse = (trivial - performance(:,5)).^2;

clear allsegs;

% the melodic distance of each note, for plotting against output curve
melDistance = Inf(size(outputmidi,1),1);
for ii = 1:size(segments,1)
    melDistance(segments{ii,2}:(segments{ii,2}+size(segments{ii,1},1)-1),1) = ...
        segments{ii,3} + zeros(size(segments{ii,1},1),1);
end
clear ii

figure
hold on
plot(mse);
plot(trivialMse);
plot(melDistance)
hold off;


% correlation between computed distance of a segment and the correlation
% value between prediction and performance of the same segment
corr_corrXmelDist = corr(cc(:,1),cc(:,3));


% correlation between the mean squared error between prediction and
% performance of a segment and the computed distance between score segment
% and reference segment.
corr_errXmelDist = corr(cc(:,2),cc(:,3));

x = cc(cc(:,3)<=2,1);
SEM = std(x)/sqrt(length(x));               % Standard Error
ts = tinv([0.025  0.975],length(x)-1);      % T-Score
PLMIN = ts*SEM;
CIcorr = mean(x) + PLMIN;                       % Confidence Interval

x = cc(cc(:,3)<=2,2);
SEM = std(x)/sqrt(length(x));               % Standard Error
ts = tinv([0.025  0.975],length(x)-1);      % T-Score
PLMIN = ts*SEM;
CImse = mean(x) + PLMIN;                       % Confidence Interval

clear x SEM ts PLMIN folder sep