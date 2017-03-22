% Filters an audio signal according to A-weighting curves for 
% (approximately) equal-loudness perception.
function filtered = Aweight(wavfile, sRate)
    wavInF = fft(wavfile);
    n = length(wavInF);
    nyq = sRate/2;
    aWeight = arrayfun(@(f) ...
        12194^2*f^4/((f^2 + 20.6^2)*(f^2 + 12194^2)*sqrt((f^2 + 107.7^2)*(f^2 + 737.9^2))) ...
        ,-nyq:2*nyq/n:nyq);
    filteredWavInF = aWeight(1:end-1)' .* fftshift(wavInF);
    filtered = ifft(ifftshift(filteredWavInF));
end