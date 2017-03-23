% This function segments a piece based on Emilios Camboroupoulos LBDM
% with two strategies for thresholding:
% 1: threshold is twice the mean boundary probability
% 2: threshold determined by iterative gradient descent and agglutination
% or further division are applied if segments are too short or too long.
% too short = fewer than 3 notes and onset interval (duration + final 
% rest) lesser than 3 * median duration.
% too long = more than 8 notes and duration greater than 8 * median duration
function seg = findPhrasesV1(midi, method)

    too_short = 4;
    too_long = 15;
    
    function s = segmentWithThreshold(midi, bounds, threshold)
        s = {};
        segStart = 1;
        segidx = 1;
        for ii = 2:(length(midi)-1)
            if bounds(ii) > threshold
                s{segidx,1} = [midi(segStart:ii-1,:) bounds(segStart:ii-1)];
                s{segidx,2} = segStart;
                segStart = ii;
                segidx = segidx + 1;
            end
        end
        s{segidx,1} = [midi(segStart:end,:) bounds(segStart:end)];
        s{segidx,2} = segStart;
    end

    function s = segmentAndRecurse(midi, bounds, threshold, limDur)
        s = segmentWithThreshold(midi, bounds, threshold);
        % check if too large segments remain and segment them
        % with local threshold values
        i = 1;
        while (i <= length(s))
            if size(s{i,1},1) > too_long
                totDur = sum(s{i,1}(:,7));
                if (totDur > too_long * medDur)
                    % too long. recurse
                    disp('recurse with segment =');
                    disp([length([s{i,1}]),s{i,2}]);
                    saux = segmentAndRecurse(s{i,1}(:,1:7),s{i,1}(:,8),bhthres(s{i,1}(2:end,8)),limDur);
                    disp(['recurse result size ',length(saux)]);
                    saux(:,2) = num2cell(cell2mat(saux(:,2)) + (s{i,2} - 1));
                    % insert into segment array
                    s = [s(1:i-1,:); saux; s(i+1:end,:)];
                    %segtail = s(i+1:end,:);
                    %s(i:(i+length(saux)-1),:) = saux;
                    %s((i+length(saux)):(i+length(saux)+length(segtail)-1),:) = segtail;
                    i = i + length(saux) - 1;
                end
            end
            i = i + 1;
        end

        % check if there are too small segments and bind them to adjacent ones
        i = 1;
        while (i < length(s))
            while size(s{i,1},1) < too_short
                ioi = s{i+1,1}(1,6) - s{i,1}(1,6);
                if ioi < too_short * medDur
                    % too small. decide where to bind
                    if s{i,1}(1,8) > s{i+1,1}(1,8)
                        if s{i+1,1}(1,8) > thrA
                            % boundaries are too strong, skip segment
                            break
                        end
                        % else, bind with next
                        s(i,:) = {[s{i,1}; s{i+1,1}],s{i,2}};
                        s(i+1,:) = [];
                    else
                        if s{i,1}(1,8) > thrA
                            % boundaries are too strong, skip segment
                            break
                        end
                        % else, bind with previous
                        s(i-1,:) = {[s{i-1,1}; s{i,1}],s{i-1,2}};
                        s(i,:) = [];
                    end
                else
                    % large enough duration
                    break
                end
                if i == length(s) break; end
            end
            i = i+1;
        end
        % if last segment is too small...
        if size(s{end,1},1) < too_short
            totDur = sum(s{end,1}(:,7));
            if totDur < too_short * medDur
                % ... and its boundary is not too strong, ...
                if s{end,1}(1,8) < thrA
                    % bind with previous.
                    s(end-1,:) = {[s{end-1,1}; s{end,1}],s{end-1,2}};
                    s(end,:) = [];
                end
            end
        end
    end

if nargin < 2, method = 'fabio'; end
% first, run LBDM on the piece
if nargin < 3, bounds = boundary(midi); end

if strcmp(method,'sergio')
    thr = 2 * mean(bounds);
    seg = segmentWithThreshold(midi, bounds, thr);
else
    % determine threshold value
    auxBounds = bounds(3:end); % ignore bounds(1)=1 and bounds(2)=0
    thrA = bhthres(auxBounds); % highest level in hierarchy
    %thrB = bhthres(auxbounds(auxBounds < thrA)); % middle level 
    %thrC = bhthres(auxbounds(auxBounds < thrB)); % low level
    medDur = median(midi(:,7));
    % segment with top level threshold
    seg = segmentAndRecurse(midi, bounds, thrA, medDur);
end
end
