% This function sets the values of velocity in the midiin array to reflect
% the loudness given by the signal loudness provided the signal vector is
% synchronized with the onset times in the midiin array.
% loudness is assumed to be in dBFS and we convert linearly so that 
% -100 dBFS (or lower) = 0 velocity (approximate softest detectable sound in
% 16-bit pcm) and 0 dBFS = 127 velocity.
function midiout = computeVelocity(midiin, loudness)

tolerance = 2; %samples

midiout = midiin;
maxval = 0;
minval = 127;
maxloud = max(loudness(:,2));
for r = 1:size(midiin,1)
    win_open = find(loudness(:,1) < midiin(r,6), 1, 'last') - tolerance;
    win_close = find(loudness(:,1) > (midiin(r,6) + midiin(r,7)),1) + tolerance;
    
    midiout(r,5) = max(0,mean(loudness(win_open:win_close,2)) - maxloud + 100)*127/100;
    maxval = max(midiout(r,5),maxval);
    minval = min(midiout(r,5),minval);
end
