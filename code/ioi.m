function ioi_vec = ioi(score_nmat, onsets, beats_per_measure)

n = size(score_nmat,1);
ioi_vec = zeros(n,1);

score_nmat(:,1) = score_nmat(:,1)/beats_per_measure;
% first_in_measure = [measure, index of 1st note, score-based onset time,
%                     performed onset time]
first_in_measure(:,2) = find(floor(score_nmat(:,1)) - [0; floor(score_nmat(1:end-1,1))] > 0);
first_in_measure(:,1) = floor(score_nmat(first_in_measure(:,2),1));
first_in_measure(:,3) = score_nmat(first_in_measure(:,2),6);
first_in_measure(:,4) = onsets(first_in_measure(:,2),1);

for i = 1:n
    measure = floor(score_nmat(i,1));
    m = find(first_in_measure(:,1) == measure); % index of first_in_measure for first note in measure of i
    m1 = min(m+1, size(first_in_measure, 1));
    m2 = min(m+2, size(first_in_measure, 1));
    m_1 = max(1, m-1);
    m_2 = max (1, m-2);
    if first_in_measure(m) == i
        % is 1st in measure
        if m == 1 || m == size(first_in_measure,1), ioi_vec(i) = NaN; % no measurement for first or last synchronization notes
        elseif m == 2 || m == size(first_in_measure,1) - 1
            % for these notes we use a single, three-bar long synchronization window from m-1 to m+1
            ioi_vec(i) = lin_interpolation(score_nmat(i,6), ...
                first_in_measure(m_1,4), first_in_measure(m_1,3), ...
                first_in_measure(m1,4), first_in_measure(m1,3));
        else
            % for these notes we calculate the onsets based on four
            % windows: [m-2,m+2],[m-1,m+1],[m-2,m+1],[m-1,m+2]
            w = lin_interpolation( ...
                ones(4,1)*score_nmat(i,6), ...
                first_in_measure([m_2; m_1; m_2; m_1], 4), ...
                first_in_measure([m_2; m_1; m_2; m_1], 3), ...
                first_in_measure([m2; m1; m1; m2], 4), ...
                first_in_measure([m2; m1; m1; m2], 3));
            % and the output is the mean of the two central values 
            % (minimum and maximum are discarded)
            ioi_vec(i) = (sum(w) - min(w) - max(w))/2;
        end
        % all other notes are synced with the first note of its measure and
        % the first note of the next measure
    else
        % except for notes in the last measure (since there is no "next
        % measure" for them) which use the synchronization notes from the
        % previous measure
        if (m == m1), m = m_1; end
        ioi_vec(i) = lin_interpolation(score_nmat(i,6), ...
            first_in_measure(m, 4), first_in_measure(m, 3), ...
            first_in_measure(m1, 4), first_in_measure(m1, 3));
    end
end

end