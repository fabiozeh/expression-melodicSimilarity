% This function sets the values of velocity in the midiin array to reflect
% the loudness in dBFS computed from the wavfile, provided it is
% synchronized with the onset times in the midiin array.
function midiout = computeNoteLoudness(midiin, wavfile, sRate)

midiout = midiin;
midiout(:,5) = meanLoudnessdBFS(midiin(:,6), midiin(:,6)+midiin(:,7), wavfile, sRate);
midiout = dbfs2vel(midiout);
end