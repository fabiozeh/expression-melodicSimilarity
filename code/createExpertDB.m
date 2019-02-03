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
%   entire piece.
%   6: z-score of mean loudness level in segment, relative to piece mean
%   7: duration of segment in seconds

% input can be a cell array of wav file names or a single folder name.
function [expertDB, params] = createExpertDB(input, applyAweighting, isTrainingSet)

if nargin < 2, applyAweighting = 0; end
if nargin < 3, isTrainingSet = 0; end

if isunix(), sep = '/'; else sep = '\'; end

% Load the miditoolbox library
addpath(genpath('miditoolbox'));

% Load the musicxml parser
addpath(genpath('util/musicxml'));

expertDB = {};
params = [];

if (~iscell(input))
    d = dir(folder);
    d = d(contains({d.name},'.wav'));
    pieceList = {d.name};
    pieceList = strcat([folder sep], pieceList);
else
    pieceList = input;
end
for piece = pieceList

    % collect score and expressive features for this piece
    [score, alignedperf] = exprFeat(piece{1}, applyAweighting, 1);
    scorefeats = score(:,8:end);
    score = score(:,1:7);
    alignedperf(isnan(alignedperf(:,9)),9) = 0; % default onset deviation = 0
    
    % segment the score into melodic phrases
    if (isTrainingSet)
        % segment according to performance timing (for better results)
        segments = findPhrases(score, 10, 0, boundary(alignedperf), 1); % segmentOnMeasures(score, 4, 1, 2);
    else
        % if it's a test set, we want the segments as they are obtained from
        % score
        segments = findPhrases(score);
    end

    % compute piece overall mean dynamics and onset-dev.
    overall_dyn = alignedperf(:,5)'*alignedperf(:,7)./sum(alignedperf(:,7));
    overall_odev = mean(alignedperf(:,9));
    
    % compute piece dynamic range and onset-dev. std.
    dynRange = max(alignedperf(:,5)) - min(alignedperf(:,5));%std(alignedperf(:,5), alignedperf(:,7));
    odevStd = std(alignedperf(:,9));
    
    % get piece tempo
    tempo = mean(alignedperf(:,8));
    
    % copy performance information into segment cell array
    segments(:,3:7) = num2cell(deal(0)); % preallocation
    for ii = 1:size(segments,1)
        % performance dynamics as midi velocity
        segdyn = alignedperf(segments{ii,2}:(segments{ii,2} + size(segments{ii,1},1) - 1),[5,6,7,1,2,9,8]);
        segments{ii,1}(:,5) = segdyn(:,1);
        % performance timing
        segments{ii,1}(:,9:12) = segdyn(:,2:5);
        % piece name
        segments{ii,3} = piece;
        % mean velocity
        segments{ii,4} = segdyn(:,1)'*segdyn(:,3)./sum(segdyn(:,3));
        % mean velocity offset from piece mean (salience)
        segments{ii,5} = segments{ii,4} - overall_dyn;
        % duration of segment in seconds
        segments{ii,7} = segments{ii,1}(end,9) + segments{ii,1}(end,10) - segments{ii,1}(1,9);
        % alpha ((mean - overall) / range);
        segments{ii,8} = segments{ii,5}./dynRange;
        % beta (segment range (std) divided by piece range)
        segments{ii,9} = (max(segdyn(:,1))-min(segdyn(:,1)))./dynRange;%std(segdyn(:,1),segdyn(:,3))./dynRange;
        % gamma (contour z-score)
        if size(segments{ii,1},1) < 2
            segments{ii,10} = 0;
        else
            segments{ii,10} = (segdyn(:,1) - segments{ii,4})./(dynRange.*segments{ii,9});
        end
        % normalized onset deviations
        segments{ii,11} = (segdyn(:,6) - overall_odev) ./ odevStd;
        % local tempo
        segments{ii,1}(:,13) = segdyn(:,7)./tempo;
    end
    
    % segment velocity z-score
    segments(:,6) = num2cell(zscore([segments{:,4}]));
    
    % include all computed data into expertDB output cell array
    expertDB = [expertDB; segments]; %#ok<AGROW>
    params = [params; overall_dyn, dynRange, overall_odev, odevStd]; %#ok<AGROW>
end

end
