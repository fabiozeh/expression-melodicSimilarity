function l = meanLoudnessdBFS(start_time, end_time, wavfile, sRate)
if ~iscolumn(start_time), start_time = start_time'; end
if ~iscolumn(end_time), end_time = end_time'; end

tolerance = 2;
indexStart = max(0, floor(start_time .* sRate) - tolerance);
indexEnd = min(size(wavfile,1), ceil(end_time .* sRate) + tolerance);

s = min(size(indexStart,1),size(indexEnd,1));
l = zeros(s,1);
for i = 1:s
    l(i) = 20*log10(rms(wavfile(indexStart(i):indexEnd(i))));
end

end