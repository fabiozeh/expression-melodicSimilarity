% Perform onset deviation estimation for a given score.
% Parameters:
% m: the desired mean deviation (seconds)
% d: the desired standard deviation (seconds)
% expertDB: the training set, created with createExpertDB.m
function output = onsetDevEstimation(scoremidi, m, d, expertDB)


%% segment the score into melodic phrases
segments = findPhrases(scoremidi);
S = size(expertDB, 1);
s = size(segments,1);

%% Compute melodic distances and auxiliary values for each segment
scores = Inf(s,S);
segments(:,3:4) = num2cell(deal(0));
for i = 1:s
    for j = 1:S
        [H, tbk] = dtwSig2(segments{i,1}, expertDB{j,1}, 0, 1, 0, 1, 'no');
        scores(i,j) = H(tbk(end,1),tbk(end,2));
    end
    x0 = segments{i,1}(1,1);
    x1 = segments{i,1}(end,1) + segments{i,1}(end,2);
    % normalized x points for output curve reading
    segments{i,3} = lin_interpolation(segments{i,1}(:,1) + 0.5*segments{i,1}(:,2), ...
        0, x0, 10, x1);
end

% compute the quadratic coefficients that approximate the contour curve for
% each expertDB segment

quadCoef = zeros(S,3);
expertDBxs = cell(S,1);
for i = 1:S
    x = expertDB{i,1};
    x0 = x(1,1);
    x1 = x(end,1)+x(end,2);
    x = x(:,1)+0.5.*x(:,2);
    expertDBxs{i} = lin_interpolation(x, 0, x0, 10, x1);
    quadCoef(i,:) = polyfit(expertDBxs{i}, ...
        expertDB{i,11},2);
end

clear H tbk x0 x1

%% Nearest neighbor estimates

nn_calc = cell(s,8);
[nn_calc{:,1}] = deal(Inf);
for i = 1:s
    for j = 1:S
        if scores(i,j) < nn_calc{i,1}
            nn_calc{i,1} = scores(i,j); % computed melodic distance score
            nn_calc{i,2} = j; % index of match(es) in expertDB
        elseif scores(i,j) == nn_calc{i,1}
            nn_calc{i,2} = [nn_calc{i,2}; j];
        end
    end
    C = quadCoef(nn_calc{i,2}(1),:); % nearest neighbor quad coefficients
    nn_calc{i,3} = C(1).*segments{i,3}.^2 + C(2).*segments{i,3} + C(3); % contour
end

clear i j matchSeg x x0 x1

% prepare output
output = segments;
output(:,3) = nn_calc(:,3);
for i = 1:s
    output{i,1}(:,6) = output{i,1}(:,6) + d .*(m + output{i,3});
    output{i,1}(:,7) = output{i,1}(:,7) - d .*(m + output{i,3});
end


end