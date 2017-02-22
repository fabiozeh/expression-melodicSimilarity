% This function segments a piece based on Emilios Camboroupoulos LBDM
% so that each segment has at most max_notes notes unless its duration
% is shorter than lim_dur.
function segments = findPhrases(midi, max_notes, lim_dur, bounds, offset)

if nargin < 2
    max_notes = 10;
    lim_dur = 0; % median(midi(:,7));
    bounds = boundary(midi);
    offset = 1;
end;

if length(midi) <= max_notes
    ind = 1; % no split
else
    totDur = sum(midi(:,7));
    if (totDur < max_notes * lim_dur)
        ind = 1;
    else
    % if there is a locally large (z-score > 2) boundary, split
    [t, ind] = max(zscore(bounds(3:end-1))); % bounds(1) is always largest
                                       % and we avoid splitting into
                                       % segments of only one note.
    if t > 2, ind = ind + 2; else ind = 1; end
    end
end

if ind == 1
    segments = {[midi, bounds], offset}; %one segment
else
    segments = [
        findPhrases(midi(1:ind-1,:), max_notes, lim_dur, bounds(1:ind-1,:), offset);
        findPhrases(midi(ind:end,:), max_notes, lim_dur, bounds(ind:end,:), offset+ind-1)
        ];
end

end
