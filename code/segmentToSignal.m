% The function takes in a midi array nmat and transforms it into an array
% representing a discrete time series of notes.
function fakeSig = segmentToSignal(nmat)

% determine the minimum multiplication factor that can represent onset 
% timing in an integer number of samples.
function k = multFactor(x)
[~, k] = rat(x);
end
f = lcm(sym(arrayfun(@(x) multFactor(x), nmat(:,1))));
% multiply onset timing by factor
nmat(:,1) = f.*nmat(:,1);
fakeSig = zeros(int32(nmat(end,1)+ceil(nmat(end,2)*f)));
sigInd = 1;
for ind = 1:(size(nmat,1)-1)
    %create samples of next note
    n = int32(nmat(ind+1,1) - nmat(ind,1));
    fakeSig(sigInd:sigInd+n) = nmat(ind,4).*ones(n);
end
n = int32(ceil(nmat(end,2)*f));
fakeSig(sigInd:sigInd+n) = nmat(end,4).*ones(n);
end
