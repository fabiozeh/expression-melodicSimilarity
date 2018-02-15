function f = midiPitchToFreq(p) 

c1 = 8.176;
factor = 1.059463094359295;

f = c1 * factor ^ p;