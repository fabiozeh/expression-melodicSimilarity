function psig = pitchSignal(nmat, sRate)

addpath(genpath('util'));

psig = zeros(ceil((nmat(end, 6) + nmat(end, 7)) * sRate), 1);
for i = 1:size(nmat,1)

    p = midiPitchToFreq(nmat(i,4));
    st = 1+round(nmat(i,6)*sRate);
    e = st+round(nmat(i,7)*sRate);
    psig(st:e) = p;
    
end
    
end