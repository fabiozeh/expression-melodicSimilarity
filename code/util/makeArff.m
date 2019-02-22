function makeArff (name, nmat)

if isunix(), sep = '/'; else, sep = '\'; end
addpath(genpath(['..' sep 'miditoolbox']));

features = nmat(:,1:19);

features(:,20:22) = circshift(features(:,[2 4 7]),1); % previous duration_beats, pitch and duration_s
features(1,20:22) = NaN; % no previous value for first note

features(:,23) = features(:,2) - features(:,21); % interval from previous note

features(:,24:27) = circshift(features(:,[2 4 7 23]),-1); % next duration_beats, pitch, duration_s and interval
features(end:24:27) = NaN; % no next value for last notes

features(:,28) = rem(features(:,4),12); % pitch, disregarding octave (0 = C, 1 = C#, ...)

% key-based pitch (e.g.: 0 for tonic, 7 for fifth)
key = mod(features(:,9)*7 + 12,12);
features(:,29) = mod(features(:,28) - key + 12, 12);

% is note a chord tone from tonic chord or relative minor?
features(:,30) = ismember(features(:,29), [0 4 7 9]);

% is note a chord tone from the V chord?
features(:,31) = ismember(features(:,29), [2 7 11]);

% chord probabilities
features(:,32:38) = chordProbabilities(features(:,29),features(:,8));

% chord with max probability (0 = tonic, 1 = supertonic, 2 = mediant...)
[~,features(:,39)] = max(features(:,32:38),[],2);
features(:,39) = features(:,39) - 1;

% note interval to likely chord root
roots = [0 2 4 5 7 9 11];
features(:,40) = mod(features(:,29) - roots(features(:,39)+1)' + 12, 12);

% is note a chord tone?
features(:,41) = (ismember(features(:,39), [0 3 4]) & ismember(features(:,40), [0 4 7])) | ...
    (ismember(features(:,39), [1 2 5]) & ismember(features(:,40), [0 3 7])) | ...
    (ismember(features(:,39), 6) & ismember(features(:,40), [0 3 6]));

%Narmour (MIDI toolbox)
features(:, 42) = narmour(nmat,'rd'); %registral direction (revised, Schellenberg 1997)
features(:, 43) = narmour(nmat,'rr'); %registral return (revised, Schellenberg 1997)
features(:, 44) = narmour(nmat,'id'); %intervallic difference
features(:, 45) = narmour(nmat,'cl'); %closure
features(:, 46) = narmour(nmat,'pr'); %proximity (revised, Schellenberg 1997)
features(:, 47) = narmour(nmat,'co'); %consonance (Krumhansl, 1995)

%features that require looping
for i = 1:size(features,1)
    features(i,48) = rem(features(i,1), features(i,10)); % onset beat (within measure)
    % metric strength
    if (features(i,48) == 0)
        features(i,49) = 3;
    elseif (features(i,48) == 2)
        features(i,49) = 2;
    elseif (features(i,48) - floor(features(i,48)) == 0)
        features(i,49) = 1;
    else
        features(i,49) = 0;
    end     
end

header = ['@relation ' name newline ...
newline ...
'@attribute duration_beats numeric' newline ...
'@attribute pitch numeric' newline ...
'@attribute duration_s numeric' newline ...
'@attribute measure numeric' newline ...
'@attribute dynamics {n,2,3,4,5,6,7,8,9,s,0}' newline ...
'@attribute beats_since_dyn numeric' newline ...
'@attribute dyn_change {n,c,d}' newline ...
'@attribute articulation {l,.,<,-}' newline ...
'@attribute score_vib {0,1}' newline ...
'@attribute ornamt {n,g,t}' newline ...
'@attribute slur {0,1}' newline ...
'@attribute slur_start {0,1}' newline ...
'@attribute prev_dur_beats numeric' newline ...
'@attribute prev_pitch numeric' newline ...
'@attribute prev_dur_s numeric' newline ...
'@attribute interval_prev numeric' newline ...
'@attribute next_dur_beats numeric' newline ...
'@attribute next_pitch numeric' newline ...
'@attribute next_dur_s numeric' newline ...
'@attribute interval_next numeric' newline ...
'@attribute pitch_in_oct {0,1,2,3,4,5,6,7,8,9,10,11}' newline ...
'@attribute pitch_in_key {0,1,2,3,4,5,6,7,8,9,10,11}' newline ...
'@attribute tonic_chord_tone {0,1}' newline ...
'@attribute dominant_chord_tone {0,1}' newline ...
'@attribute prob_chord_I numeric' newline ...
'@attribute prob_chord_ii numeric' newline ...
'@attribute prob_chord_iii numeric' newline ...
'@attribute prob_chord_IV numeric' newline ...
'@attribute prob_chord_V numeric' newline ...
'@attribute prob_chord_vi numeric' newline ...
'@attribute prob_chord_vii numeric' newline ...
'@attribute chord {0,1,2,3,4,5,6}' newline ...
'@attribute pitch_in_chord {0,1,2,3,4,5,6,7,8,9,10,11}' newline ...
'@attribute current_chord_tone {0,1}' newline ...
'@attribute nar_rd numeric' newline ...
'@attribute nar_rr numeric' newline ...
'@attribute nar_id numeric' newline ...
'@attribute nar_cl numeric' newline ...
'@attribute nar_pr numeric' newline ...
'@attribute nar_co numeric' newline ...
'@attribute beat_in_measure numeric' newline ...
'@attribute metric_strength {0,1,2,3}' newline ...
newline ...
'@data' newline];

F = fopen(name,'wt');

fprintf(F,'%s',header);

for i = 1:size(nmat,1)
    fprintf(F,['%.6f,%d,%.6f,%d,%c,%d,%c,%c,%d,%c,%d,%d,' ...
        '%.6f,%d,%.6f,%d,%.6f,%d,%.6f,%d,%d,%d,%d,%d,' ...
        '%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%d,%d,%d,' ...
        '%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.6f,%d\n'],features(i,[2 4 7 8 12:49]));
end

fclose(F);

end

function probs = chordProbabilities(notes, measure)
S = size(notes,1);
probs(S,7) = 0; %preallocation
for i=1:S
    window = notes(ismember(measure, measure(i)));
    s = size(window,1);
    probs(i,1) = sum(ismember(window, [0 4 7]))/s;
    probs(i,2) = sum(ismember(window, [2 5 9]))/s;
    probs(i,3) = sum(ismember(window, [4 7 11]))/s;
    probs(i,4) = sum(ismember(window, [5 9 0]))/s;
    probs(i,5) = sum(ismember(window, [7 11 2]))/s;
    probs(i,6) = sum(ismember(window, [9 0 4]))/s;
    probs(i,7) = sum(ismember(window, [11 2 5]))/s;
end
end
