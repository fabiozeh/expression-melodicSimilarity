% Perform local tempo estimation for a given score.
% Parameters:
% expertDB: the training set, created with createExpertDB.m
function output = localTempoEstimation(score,expertDB,tempo)

if nargin < 3
    if ~iscell(score)
        tempo = gettempo(score);
    else
        tempo = gettempo(vertcat(score{:,1}));
    end
end

%% segment the score into melodic phrases
if ~iscell(score)
    segments = findPhrases(score);
else
    segments = score;
end
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
    curve = expertDB{i,1}(:,13);
    x = expertDB{i,1};
    x0 = x(1,1);
    x1 = x(end,1)+x(end,2);
    x = x(:,1)+0.5.*x(:,2);
    expertDBxs{i} = lin_interpolation(x, 0, x0, 10, x1);
    %exclude outliers
    expertDBxs{i} = expertDBxs{i}(curve > 0.5 & curve < 2.0);
    curve = curve(curve > 0.5 & curve < 2.0);
    % find best quadratic coeficients
    quadCoef(i,:) = polyfit(expertDBxs{i}, ...
        10.*log10(curve),2);
end

clear H tbk x0 x1

%% Nearest neighbor estimates

nn_calc = cell(s,3);
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
k=2;
w = 1./(scores+1e-6);
w2 = zeros(s,S);
for i = 1:k
    [~, ind] = max(w,[],2);
    w(sub2ind([s,S],1:s,ind')) = 0;
    w2(sub2ind([s,S],1:s,ind')) = 1;
end
        
wCoef = w2 * quadCoef ./ (sum(w2,2)*ones(1,3));
for i = 1:s
    knn{i} = wCoef(i,1).*segments{i,3}.^2 + ...
        wCoef(i,2).*segments{i,3} + wCoef(i,3);
end
     
for i = 1:s
    for j = 1:size(segments{i,1},1)
        old = segments{i,1}(j,7);
        segments{i,1}(j,7) = segments{i,1}(j,2)/(tempo.*10^(nn_calc{i,3}(j)./10))*60;
        segments{i,1}(j+1:end,6) = segments{i,1}(j+1:end,6) - old + segments{i,1}(j,7);
        for ii = (i+1):s
            segments{ii,1}(1:end,6) = segments{ii,1}(1:end,6) - old + segments{i,1}(j,7);
        end
    end
end

output(:,1) = segments(:,1);
output(:,2) = nn_calc(:,1); 
output(:,3) = nn_calc(:,2);
output(:,4) = nn_calc(:,3);
for i = 1:size(output,1)
    output{i,2} = output{i,2}./size(output{i,1},1);
end

end