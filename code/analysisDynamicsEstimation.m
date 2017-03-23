%% analysis of dynamics estimation %%

%% create expert database
createExpertDB

if isunix(), sep = '/'; else sep = '\'; end

folder = ['..' sep 'data' sep 'caprice' sep 'input'];

%% Load the required libraries and files for score and performance

addpath(genpath('miditoolbox'));

% load score and calculate expressive features from performance
[score, performance] = exprFeat(folder, 1, 0);

%% estimate dynamics for target score

segments = dynamicsEstimation(score, expertDB);

% reconstruct predicted performance from segments
outputmidi = cat(1, segments{:,1});

%% test quality of output according to original performance

% TODO rewrite analysis from main
