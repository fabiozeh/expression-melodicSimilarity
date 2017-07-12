function segments = segmentOnMeasures(midiin, beatsPerMeasure, anaCruziBeats, measuresPerSegment)

current = -measuresPerSegment;
idx = 0;
segments = {};
for i = 1:size(midiin,1)
    measure = floor((midiin(i,1)+beatsPerMeasure-anaCruziBeats)/beatsPerMeasure);
    if measure - current >= measuresPerSegment
        current = measure;
        idx = idx + 1;
        segments{idx,1} = midiin(i,:); %#ok<*AGROW>
        segments{idx,2} = i;
    else
        segments{idx,1} = [segments{idx,1}; midiin(i,:)];
    end
end

end