% Extracts the expressive features from the piece in inputFolder. Depends
% on miditoolbox and musicxml parser (already in path).
function [scoremidi, alignedperf, wavfile, sRate] = exprFeat(inputFolder, detectOnsets, applyAweighting, useVel)

    if nargin < 4, useVel = 1; end

    if isunix(), sep = '/'; else sep = '\'; end

    if detectOnsets
        % Automatic note onset detection with pYin vamp plugin (running script)
        pYinNotes(inputFolder, 'performance.wav');
        perfmidi = readmidi([inputFolder sep 'performance.mid']);
    else
        perfmidi = readmidi([inputFolder sep 'perfAlignment.mid']);
    end
    
    if exist([inputFolder sep 'score.xml'], 'file')
        scoremidi = parseMusicXML([inputFolder sep 'score.xml']);
        scoremidi = scoremidi(scoremidi(:,4) ~= 0,:); % delete rests
    elseif exist([inputFolder sep 'score.mid'], 'file')
        scoremidi = readmidi([inputFolder sep 'score.mid']);
    else
        scoremidi = perfmidi;
    end

    %% Calculate dynamics of expert performances

    % load the performance audio signal
    [wavfile, sRate] = audioread([inputFolder sep 'performance.wav']);

    % filter it for more realistic loudness calculation from an auditory
    % perspective
    if applyAweighting
        wavfile = Aweight(wavfile, sRate);
    end

    % compute average energy around each performed note
    perfmidi = computeNoteLoudness(perfmidi, wavfile, sRate, useVel);

    % align the performance with the score
    [scoremidi, alignedperf] = perfAlign(scoremidi, perfmidi);
    % decide an artificial note for deletions
    % TODO
    
    % compute local tempo
    alignedperf(:,8) = scoremidi(:,2).*60./alignedperf(:,7);
    
    % compute timing deviation
    alignedperf(:,9) = timing(alignedperf, scoremidi);
    
end