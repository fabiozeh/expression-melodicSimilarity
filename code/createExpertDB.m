% Generates cell array 'expertDB' with melody segments and expressive 
% features from a set of input files given by setting variable
% folderList appropriately.

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

clear ii segments piece score alignedperf folder sep folderList applyAweighting
