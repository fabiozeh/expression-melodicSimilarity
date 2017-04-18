% Generates cell array 'expertDB' with melody segments and expressive 
% features from a set of input files given by setting variable
% folderList appropriately.
% Structure of expertDB:
%   Each row has information on one melodic segment of a piece in the
%   database. Columns contain the following data:
%   1: array of midi notes (typical miditoolbox nmat) of the segment.
%   Midi timing information is relative to the score but onset and duration
%   from the reference performance (in seconds) are available in extra
%   columns, respectively 9 and 10.
%   2: Index of first segment note in the original piece. E.g.:
%   expertDB{i,2} = 20 means segment "i" contains a sequence of notes
%   starting from note 20 of the piece from which it originates.
%   3: name of the piece this segment was drawn from, as given by the
%   folder names input in folderList.
%   4: mean loudness level in segment.
%   5: difference between mean loudness level in this segment and in the
%   previous segment of its originating piece.

% column 1 = folder name. Column 2 = 1 to do automatic onset detection
folderList = {'telmi_eulalie', 1; 'bachManual', 0};
applyAweighting = 0;

if isunix(), sep = '/'; else sep = '\'; end

% Load the miditoolbox library
addpath(genpath('miditoolbox'));

expertDB = {};
for piece = folderList'

    folder = ['..' sep 'data' sep piece{1} sep 'input'];
    [score, alignedperf] = exprFeat(folder, piece{2}, applyAweighting);
    
    % segment the score into melodic phrases
    segments = findPhrases(score);

    % copy performance information into segment cell array
    for ii = 1:size(segments,1)
        segments{ii,1}(:,5) = alignedperf(segments{ii,2}:(segments{ii,2} + size(segments{ii,1},1) - 1),5);
        segments{ii,1}(:,9) = alignedperf(segments{ii,2}:(segments{ii,2} + size(segments{ii,1},1) - 1),6);
        segments{ii,1}(:,10) = alignedperf(segments{ii,2}:(segments{ii,2} + size(segments{ii,1},1) - 1),7);
        segments{ii,3} = piece{1};
    end
    
    expertDB = [expertDB; segments];
    
end

% Compute useful information about each segment

% segment mean loudness level
expertDB(:,4) = num2cell(cellfun(@(x) mean(x(:,5)), expertDB(:,1)));
% offset from previous segment
expertDB(:,5) = num2cell([0; diff(vertcat(expertDB{:,4}))]);
[expertDB{[expertDB{:,2}] == 1, 5}] = deal(0);

clear ii segments piece score alignedperf folder sep folderList applyAweighting
