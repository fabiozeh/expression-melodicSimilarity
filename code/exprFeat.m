function [scoremidi, alignedperf, wavfile, sRate] = exprFeat(inputFolder, detectOnsets, applyAweighting)
    
    if isunix(), sep = '/'; else sep = '\'; end

    if detectOnsets
        % Automatic note onset detection with pYin vamp plugin (running script)
        pYinNotes(inputFolder, 'performance.wav');
        perfmidi = readmidi([inputFolder sep 'performance.mid']);
    else
        perfmidi = readmidi([inputFolder sep 'perfAlignment.mid']);
    end

    if exist([inputFolder sep 'score.mid'], 'file')
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
    perfmidi = computeNoteLoudness(perfmidi, wavfile, sRate);

    % align the performance with the score
    alignedperf = perfAlign(scoremidi, perfmidi);
    % decide an artificial note for deletions
    % TODO
end