function [aligneds, alignedp] = perfAlign(scoremidi, perfmidi)
    [~, alignment] = bestalign(scoremidi(:,4),perfmidi(:,4));
    ix1 = 0; 
    ix2 = 0;
    ixS = 1; %score midi index
    ixP = 1; %perf midi index
    ix3 = 1; %aligned arrays index
    alignedp = NaN(size(scoremidi));
    aligneds = NaN(size(scoremidi));
    for ixal = 1:size(alignment,1)
        if alignment(ixal,1) - ix1 == 0
            % insertion
            ixP = ixP + 1;
        elseif alignment(ixal,2) - ix2 == 0
            % deletion
            alignedp(ix3,:) = []; % delete one NA element
            aligneds(ix3,:) = []; % delete one NA element
            ixS = ixS + 1;
        else
            % match or error
            alignedp(ix3,:) = perfmidi(ixP,:);
            aligneds(ix3,:) = scoremidi(ixS,:);
            ix3 = ix3 + 1;
            ixS = ixS + 1;
            ixP = ixP + 1;
        end
        ix1 = alignment(ixal,1);
        ix2 = alignment(ixal,2);
    end 
end
