% Extracts the expressive features from the piece in inputFolder. Depends
% on miditoolbox and musicxml parser (already in path).
function [scoremidi, alignedperf, wavfile, sRate] = exprFeat(filename, applyAweighting, useVel)

    if nargin < 3, useVel = 1; end

    % remove filename extension
    if (filename(end-3) == '.')
        filename = filename(1:end-3);
    end
    
    if ~exist([filename 'mid'], 'file')
        if ~exist([filename 'midi'], 'file')
            % Automatic note onset detection with pYin vamp plugin (running script)
            % pYinNotes('.', filename);
            % perfmidi = readmidi([filename 'mid']);
            error('No performance alignment found.');
        else
            perfmidi = readmidi([filename 'midi']);
        end
    else
        perfmidi = readmidi([filename 'mid']);
    end
    
    
    if exist([filename 'xml'], 'file')
        scoremidi = parseMusicXML([filename 'xml']);
        scoremidi = scoremidi(scoremidi(:,4) ~= 0,:); % delete rests
    else
        scoremidi = perfmidi;
    end

    %% Calculate dynamics of expert performances

    % load the performance audio signal
    [wavfile, sRate] = audioread([filename 'wav']);

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
    
    % tempo adjustments
    %scoremidi(:,1) = scoremidi(:,1) + alignedperf(1,1) - scoremidi(1,1);
    %scoremidi(:,6) = scoremidi(:,6) + alignedperf(1,6) - scoremidi(1,6);
    alignedperf(:,1:2) = scoremidi(:,1:2);
    
    % compute local tempo
    alignedperf(:,8) = alignedperf(:,2).*60./alignedperf(:,7);
    
    % adjust score tempo according to performance
    scoremidi = settempo(scoremidi, gettempo(alignedperf));
    
    % compute timing deviation (if score came from xml)
    if (size(scoremidi,2) > 7)
        alignedperf(:,9) = timing(alignedperf, scoremidi);
    end
    
end