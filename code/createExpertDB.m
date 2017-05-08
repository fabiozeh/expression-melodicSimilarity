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
%   6: z-score of mean loudness level in segment, relative to piece mean

% column 1 = folder name. Column 2 = 1 to do automatic onset detection
folderList = {'meditacion', 0; 'bachManual', 0; 'beethoven4_3', 0; 'borodin2_1', 0; 'haydn', 0};
applyAweighting = 1;

if isunix(), sep = '/'; else sep = '\'; end

% Load the miditoolbox library
addpath(genpath('miditoolbox'));

expertDB = {};
for piece = folderList'

    folder = ['..' sep 'data' sep piece{1} sep 'input'];
    
    % collect score and expressive features for this piece
    [score, alignedperf] = exprFeat(folder, piece{2}, applyAweighting);
    
    % segment the score into melodic phrases
    segments = findPhrases(score);

    % copy performance information into segment cell array
    segments(:,3:5) = num2cell(deal(0)); % preallocation
    for ii = 1:size(segments,1)
        % performance dynamics as midi velocity
        segments{ii,1}(:,5) = alignedperf(segments{ii,2}:(segments{ii,2} + size(segments{ii,1},1) - 1),5);
        % performance timing
        segments{ii,1}(:,9) = alignedperf(segments{ii,2}:(segments{ii,2} + size(segments{ii,1},1) - 1),6);
        segments{ii,1}(:,10) = alignedperf(segments{ii,2}:(segments{ii,2} + size(segments{ii,1},1) - 1),7);
        % piece name
        segments{ii,3} = piece{1};
        % mean velocity
        segments{ii,4} = mean(segments{ii,1}(:,5));
        % mean velocity offset from previous segment
        if ii > 1
            segments{ii,5} = segments{ii,4} - segments{ii-1,4};
        end
    end
    
    % segment velocity z-score
    segments(:,6) = num2cell(zscore([segments{:,4}]));
    
    % include all computed data into expertDB output cell array
    expertDB = [expertDB; segments];
    
end

clear ii segments piece score alignedperf folder sep folderList applyAweighting
