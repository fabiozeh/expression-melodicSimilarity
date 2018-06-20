function t = timing(perfmat, scoremat)

s = size(perfmat, 1);
sync0 = NaN;
sync1 = NaN;
t = NaN(s,1);
firstInMeasure = [1; diff(scoremat(:,8))];
measures = scoremat(logical(firstInMeasure),8);
onsets1stInM_s = scoremat(logical(firstInMeasure),6);
onsets1stInM_p = perfmat(logical(firstInMeasure),6);
for i = 1:s
    m = find(measures == scoremat(i,8));
    m1 = min([m+1, size(measures,1)]);
    m2 = min([m+2, size(measures,1)]);
    m_1 = max([m-1, 1]);
    m_2 = max([m-2, 1]);
    if (firstInMeasure(i))
        if (m == 1 || m == size(measures,1))
            continue; % no measurement for first or last synchronization notes
        elseif (m == 2 || m == size(measures,1)-1)
            % for these notes we use a single, three-bar long synchronization window from m-1 to m+1
            t(i) = lin_interpolation(scoremat(i,6), ...
                onsets1stInM_p(m_1), onsets1stInM_s(m_1), ...
                onsets1stInM_p(m1), onsets1stInM_s(m1));
        else
            % for these notes we calculate the onsets based on four windows
            % from measure m-2 to m+2
            w5 = lin_interpolation(scoremat(i,6), ...
                onsets1stInM_p(m_2), onsets1stInM_s(m_2), ...
                onsets1stInM_p(m2), onsets1stInM_s(m2));
            % from measure m-1 to m+1
            w3 = lin_interpolation(scoremat(i,6), ...
                onsets1stInM_p(m_1), onsets1stInM_s(m_1), ...
                onsets1stInM_p(m1), onsets1stInM_s(m1));
            % from measure m-2 to m+1
            w4l = lin_interpolation(scoremat(i,6), ...
                onsets1stInM_p(m_2), onsets1stInM_s(m_2), ...
                onsets1stInM_p(m1), onsets1stInM_s(m1));
            % and from measure m-1 to m+2
            w4r = lin_interpolation(scoremat(i,6), ...
                onsets1stInM_p(m_1), onsets1stInM_s(m_1), ...
                onsets1stInM_p(m2), onsets1stInM_s(m2));
            % and the output is the average of the two central values
            % (minimum and maximum are discarded)
            t(i) = (w3 + w5 + w4r + w4l - min([w3,w5,w4r,w4l]) - max([w3,w5,w4r,w4l])) ./ 2.0;
        end
    else
        % all other notes are synced with the first note of its measure and
        % the first note of the next measure
        if (m == m1)
            % except for notes in the last measure (since there is no "next measure" for them)
            % which use the synchronization notes from the previous measure
            m = m - 1;
        end
        t(i) = lin_interpolation(scoremat(i,6), ...
            onsets1stInM_p(m), onsets1stInM_s(m), ...
            onsets1stInM_p(m1), onsets1stInM_s(m1));
    end
end

t = perfmat(:,6) - t; % ideal timing - actual timing = deviation

end