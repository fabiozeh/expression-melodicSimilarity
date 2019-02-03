% Extracts the expressive features from the piece in inputFolder. Depends
% on miditoolbox and musicxml parser (already in path).
function [scoremidi, alignedperf, wavfile, sRate] = exprFeat(filename, applyAweighting, useVel)

    if nargin < 3, useVel = 1; end

    if ~exist([filename(1:end-3) 'mid'], 'file')
        % Automatic note onset detection with pYin vamp plugin (running script)
        pYinNotes('.', filename);
    end
    
    perfmidi = readmidi([filename(1:end-3) 'mid']);
    
    if exist([filename(1:end-3) 'xml'], 'file')
        scoremidi = parseMusicXML([filename(1:end-3) 'xml']);
        scoremidi = scoremidi(scoremidi(:,4) ~= 0,:); % delete rests
        % grace note inclusion
%         aux = scoremidi;
%         ptr = 1;
%         for l = 1:size(aux,1);
%            if (aux(l,17) == 'g') % grace note
%                scoremidi = [scoremidi(1:ptr-1,:); scoremidi(ptr,:); scoremidi(ptr:end,:)];
%                dur = 0;
%                if (scoremidi(ptr,2) >= 0.5)
%                    dur = 0.25;
%                else
%                    dur = scoremidi(ptr,2)/2;
%                end
%                    lcltmp = scoremidi(ptr,7)/scoremidi(ptr,2);
%                    scoremidi(ptr-1,[2,7]) = [scoremidi(ptr-1,2)-dur, scoremidi(ptr-1,7)-dur*lcltmp];
%                    scoremidi(ptr,[1,2,4,6,7]) = [scoremidi(ptr,1)-dur, dur, scoremidi(ptr,4)+2, ...
%                        scoremidi(ptr,6) - dur*lcltmp, dur*lcltmp];
%                ptr = ptr + 2;
%            else
%                ptr = ptr + 1;
%            end
%         end
    else
        scoremidi = perfmidi;
    end

    %% Calculate dynamics of expert performances

    % load the performance audio signal
    [wavfile, sRate] = audioread(filename);

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
    
    % compute timing deviation
    alignedperf(:,9) = timing(alignedperf, scoremidi);
    
end