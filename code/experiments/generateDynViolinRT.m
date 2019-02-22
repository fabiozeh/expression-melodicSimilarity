addpath(genpath('miditoolbox'));
nmat = readmidi('../../MML2017/CumbancheroExp.mid');
dynsig = segmentToSignal(nmat);
for i = 1:size(dynsig,2)
    dynsigb(3*i-2) = dynsig(i);
    dynsigb(3*i-1) = dynsig(i);
    dynsigb(3*i) = dynsig(i);
end
dynsig = dynsigb/(127*4);
clear dynsigb i
