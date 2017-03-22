function r = windowedRms(wavfile, window, step)
% do rms for every sRate/100 elements (441 in 44.1kHz wav file)
ptr = 1;
r = zeros(1+floor(length(wavfile)/window),1);
for ii = 1:(length(r)-1)
    r(ii,1) = rms(wavfile(ptr:(ptr+window-1)));
    ptr = ptr + step;
end
r(ii+1,1) = rms(wavfile(ptr:end));
end