% Generates cell array 'db' with note-level score and expressive 
% features from a set of input files given by setting variable
% input appropriately.
% Structure of db:
%   Each row has information on one note in the
%   database. Columns contain the following data:
%   1-42: see noteFeatures function.
%   43: phrase number, one-based, as segmented by our method
%   44: velocity number computed from performance
%   45: onset time in seconds, as measured in performance
%   46: duration in seconds, as measured in performance
%   47: local tempo in bpm
%   48: timing deviation, in seconds
%   49: piece name
% input can be a cell array of wav file names or a single folder name.
function [db, arff_info] = expressionDataset(input)

if isunix(), sep = '/'; else sep = '\'; end

% Load the miditoolbox library
addpath(genpath('miditoolbox'));
% TODO: change to other library


% Load the musicxml parser
addpath(genpath('util/musicxml'));

db = {};
params = [];

if (~iscell(input))
    d = dir(input);
    d = d(contains({d.name},'.wav'));
    pieceList = {d.name};
    pieceList = strcat([input sep], pieceList);
else
    pieceList = input;
end
n = 1;
for piece = pieceList

    % collect score and expressive features for this piece
    [score, alignedperf] = exprFeat(piece{1}, 0, 1);
    alignedperf(isnan(alignedperf(:,9)),9) = 0; % default onset deviation = 0
    
    % segment the score into melodic phrases
    % TODO - start segmentation with score information
    segments = findPhrases(score);

    % compute note-level score features
    [dbi, arff_info] = noteFeatures(score);
    
    % aggregate segmentation information as a feature
    coln = size(dbi,2) + 1;
    j = 1;
    for i = size(segments,1)
        s = size(segments{i,1},1);
        dbi(j:j+s-1,coln) = i; % phrase number
        j = j + s - 1;
    end
    
    % copy performance information into database
    coln = coln + 1;
    dbi(:,coln:coln + 4) = alignedperf(:,5:9);
    
    % include piece index in the end
    dbi(:,coln+5) = n;
    n = n+1;
    
    % add data to output array
    db = [db; dbi]; %#ok<AGROW>
end

arff_info(end+1:end+7,:) = {
    'phrase_num', 'numeric', '%d'; ...
    'velocity', 'numeric', '%d'; ...
    'perf_onset_s', 'numeric', '%.6f'; ...
    'perf_dur_s', 'numeric', '%.6f'; ...
    'local_tempo', 'numeric', '%.6f'; ...
    'timing_dev', 'numeric', '%.6f'; ...
    'piece_idx', 'numeric', '%d'; ...
    };

end
