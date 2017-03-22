function aligned = perfAlign(scoremidi, perfmidi, fillDeletions)
    score = bestalign(scoremidi(:,4),perfmidi(:,4));
    allpaths = score{2};
    isDone = 0;
    ix1 = size(allpaths,1);
    ix2 = size(allpaths,2);
    path = {};
    ixp = 1;
    while ~isDone
        path{1,ixp} = ix1-1;
        path{2,ixp} = ix2-1;
        ix1 = allpaths(ix1, ix2, 1);
        ix2 = allpaths(ix1, ix2, 2);
        if ix1 == 1 && ix2 == 1, isDone = 1; end
        ixp = ixp + 1;
    end
    alignment = cell2mat(path(:,end:-1:1));
    ix1 = 0;
    ix2 = 0;
    ixp = 1;
    outcomes = zeros(size(alignment,2));
    for ixp = 1:size(alignment,2)
        if alignment(1,ixp) - ix1 == 0
            % insertion
            outcomes(ixp) = 'i';
        elseif alignment(2,ixp) - ix2 == 0
            % deletion
            outcomes(ixp) = 'd';
        else
            % match or error
            outcomes(ixp) = 'm';
        end
        ix1 = alignment(1,ixp);
        ix2 = alignment(2,ixp);
    end
    
    aligned = NaN(size(scoremidi));
    ixp = 1;
    ix2 = 1;
    ix1 = 1;
    isDone = 0;
    while ~isDone
        if outcomes(ixp) == 'm'
            aligned(ix1,:) = perfmidi(ix2,:);
            ix2 = ix2 + 1;
            ix1 = ix1 + 1;
            ixp = ixp + 1;
        elseif outcomes(ixp) == 'd'
            if fillDeletions
                % make up a note at
            end
            ix1 = ix1 + 1;
            ixp = ixp + 1;
        elseif outcomes(ixp) == 'i'
            ixp = ixp + 1;
            ix2 = ix2 + 1;
        end
        isDone = ix1 > size(scoremidi,1);
        isDone = isDone || ix2 > size(perfmidi,1);
    end
end
