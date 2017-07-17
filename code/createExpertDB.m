% Generates cell array 'expertDB' with melody segments and expressive 
% features from a set of input files given by variable folderList, which
% should be a cell array indicating one piece per row, with columns:
% 1: piece input folder name as string (assumed to be under ../data/)
% 2: set 1 to detect onsets automatically or 0 to use an alignment file
% provided
% Structure of output expertDB:
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
%   8: "alpha" component of dynamics model (relative phrase salience)
%   9: "beta" component of dynamics model (relative phrase range)
%  10: "Gamma" component of dynamics model (relative phrase contour)
% Structure of output params:
%   Each row contains measured parameters for one piece (a row of
%   folderList). Columns are:
%   1: mean overall dynamics level
%   2: dynamic range measure (weighted standard deviation of note dynamics)

% column 1 = folder name. Column 2 = 1 to do automatic onset detection
function [expertDB, params] = createExpertDB(folderList, applyAweighting)

if nargin < 2, applyAweighting = 0; end

if isunix(), sep = '/'; else sep = '\'; end

% Load the miditoolbox library
addpath(genpath('miditoolbox'));

expertDB = {};
params = [];
for piece = folderList'

    folder = ['..' sep 'data' sep piece{1} sep 'input'];
    
    % collect score and expressive features for this piece
    [score, alignedperf] = exprFeat(folder, piece{2}, applyAweighting);
    
    % segment the score into melodic phrases
    segments = findPhrases(score); % segmentOnMeasures(score, 4, 1, 2);

    % compute piece overall mean dynamics
    overall = alignedperf(:,5)'*alignedperf(:,7)./sum(alignedperf(:,7));
    
    % compute piece dynamic range
    dynRange = std(alignedperf(:,5), alignedperf(:,7));
    
    % copy performance information into segment cell array
    segments(:,3:7) = num2cell(deal(0)); % preallocation
    for ii = 1:size(segments,1)
        % performance dynamics as midi velocity
        segdyn = alignedperf(segments{ii,2}:(segments{ii,2} + size(segments{ii,1},1) - 1),5:7);
        segments{ii,1}(:,5) = segdyn(:,1);
        % performance timing
        segments{ii,1}(:,9) = segdyn(:,2);
        segments{ii,1}(:,10) = segdyn(:,3);
        % piece name
        segments{ii,3} = piece{1};
        % mean velocity
        segments{ii,4} = segdyn(:,1)'*segdyn(:,3)./sum(segdyn(:,3));
        % mean velocity offset from piece mean (salience)
        segments{ii,5} = segments{ii,4} - overall;
        % duration of segment in seconds
        segments{ii,7} = segments{ii,1}(end,9) + segments{ii,1}(end,10) - segments{ii,1}(1,9);
        % alpha ((mean - overall) / range);
        segments{ii,8} = segments{ii,5}./dynRange;
        % beta (segment range (std) divided by piece range)
        segments{ii,9} = std(segdyn(:,1),segdyn(:,3))./dynRange;
        % gamma (contour z-score)
        if size(segments{ii,1},1) < 2
            segments{ii,10} = 0;
        else
            segments{ii,10} = (segdyn(:,1) - segments{ii,4})./(dynRange.*segments{ii,9});
        end
    end
    
    % segment velocity z-score
    segments(:,6) = num2cell(zscore([segments{:,4}]));
    
    % include all computed data into expertDB output cell array
    expertDB = [expertDB; segments]; %#ok<AGROW>
    params = [params; overall, dynRange]; %#ok<AGROW>
end

end
