% Perform dynamics estimation for a given score %
function segments = dynamicsEstimation(scoremidi, expertDB)

meanVel = {expertDB(expertDB{:,2} == 1, 3), -1};
for i = 1:size(meanVel,1)
    meanVel{i,2} = mean(cellfun(@(x) mean(x(:,5)), expertDB(expertDB{:,3}==meanVel{i,1},1)));
end

%% segment the score into melodic phrases
segments = findPhrases(scoremidi);

%% Compute dynamics estimation for each segment
segments{:,3} = Inf(size(segments, 1),1);
for jj = 1:size(segments,1)
    % find most similar segment(s)
    for ii = 1:size(expertDB,1)
        [H, tbk] = dtwSig2(segments{jj,1}, expertDB{ii,1}, 0, 1, 0, 1, 'no');
        score = H(tbk(end,1),tbk(end,2));
        if score < segments{jj,3}
            segments{jj,3} = score;
            segments{jj,4} = ii;
        elseif score == segments{jj,3}
            segments{jj,4} = [segments{jj,4}; ii];
        end
    end
    matchSeg = expertDB(segments{jj,4}(1),:);
    segments{jj,5} = mean(matchSeg{1}(:,5)); %segment mean level
    segments{jj,1}(:,5) = scaleInterpolate(size(segments{jj,1},1), matchSeg{1}(:,5) - segments{jj,5});
    %segments{jj,6} --> offset from prev. segment
    if matchSeg{2} > 1
        segments{jj,6} = segments{jj,5} - mean(expertDB{segments{jj,4}(1)-1,1}(:,5))
    else
        segments{jj,6} = 0;
    end
    segments{jj,7} = meanVel{[meanVel{:,1}]==matchSeg{3},2}; % match piece mean level
end
% should interpolation be done in time domain instead of note domain?

clear ii jj H tbk

%% 

% TODO
% more stable algorithm:
% if segment level is > 1 sd from mean, copy neighbor gradient
% else copy deviation from mean

% adjust mean level of segments
% the strategy is as follows:
% - two parameters are saved for each reference segment: the mean level of
% the piece it belongs to (pieceVel = segments{:,7}) and the offset of the segment mean relative to
% the mean level of previous segment (os = segments{:,6}).
% - after determining all predictions, the final piece mean level is adjusted to
% be the mean(pieceVel) taking all segments into account and the offset of each
% segment is corrected to os = os - mean(os).
meanOS = mean([segments{:,6}]);
predVel = mean([segments{:,7}]);
offset = predVel + segments{1,6} - meanOS;
segments{1,1}(:,5) = segments{1,1}(:,5) + offset;
for ii = 2:size(segments,1)
    offset = segments{ii,6} + offset - meanOS;
    segments{ii,1}(:,5) = segments{ii,1}(:,5) + offset;
end

clear ii ind perfStartInd perfEndInd meanOS mL meanVel offset

end