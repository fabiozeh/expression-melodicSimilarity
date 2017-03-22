% This function sets the values of velocity in the midiin array to reflect
% the loudness given by the signal rms provided the signal vector is
% synchronized with the onset times in the midiin array.
function midiout = computeVelocity(midiin, rms)

tolerance = 2; %samples

midiout = midiin;
maxval = 0;
minval = 127;
for r = 1:size(midiin,1)
    win_open = find(rms(:,1) < midiin(r,6), 1, 'last') - tolerance;
    win_close = find(rms(:,1) > (midiin(r,6) + midiin(r,7)),1) + tolerance;
    
    midiout(r,5) = mean(rms(win_open:win_close,2)) * 127;
    maxval = max(midiout(r,5),maxval);
    minval = min(midiout(r,5),minval);
end

% center around 63.5
midiout(:,5) = midiout(:,5) + 63.5 - (maxval - minval)/2;
