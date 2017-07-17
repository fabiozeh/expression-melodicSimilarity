function segments = segmentOnMeasures(midiin, beatsPerMeasure, anaCruziBeats, measuresPerSegment)

current = 0;
startOffset = 0;
idx = 1;
segments{1,1} = [];
segments{1,2} = 1;
for i = 1:size(midiin,1)
    measure = floor((midiin(i,1)+beatsPerMeasure-anaCruziBeats)/beatsPerMeasure);
    if measure - current >= measuresPerSegment - startOffset
        startOffset = 0;
        current = measure;
        idx = idx + 1;
        segments{idx,1} = midiin(i,:); %#ok<*AGROW>
        segments{idx,2} = i;
    else
        segments{idx,1} = [segments{idx,1}; midiin(i,:)];
    end
end

end