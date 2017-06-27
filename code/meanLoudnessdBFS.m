function l = meanLoudnessdBFS(start_time, end_time, wavfile, sRate)
if ~iscolumn(start_time), start_time = start_time'; end
if ~iscolumn(end_time), end_time = end_time'; end

indexStart = start_time .* sRate;
indexEnd = end_time .* sRate;

s = min(size(indexStart,1),size(indexEnd,1));
l = zeros(s,1);
for i = 1:s
    l(i) = 20*log(rms(wavfile(indexStart(i):indexEnd(i))));
end

end