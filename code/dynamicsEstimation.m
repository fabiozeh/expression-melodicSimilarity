% Perform dynamics estimation for a given score
function [output, outputB, outputC, outputD] = dynamicsEstimation(scoremidi, expertDB, method, k, c)

allmethods = {'w-knn', 'knn', 'nn-quad', 'nn'};
m = zeros(1,size(allmethods,1));
for i = 1:size(method, 1)
	m = m | strcmp(method(i,:), allmethods);
end
if sum(m) < 1, m(1) = 1; end
if strcmp(method, 'all'), m(:) = 1; end
method = m;
idxwknn = 1;
idxknn = 2;
idxnnq = 3;
idxnn = 4;
clear m allmethods

if (method(idxknn) || method(idxwknn)) && nargin < 4, k = 7; c = 1e-9; end
if (method(idxknn) || method(idxwknn)) && nargin < 5, c = 1e-9; end

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
    % durations of segments in secs
    segments{i,3} = segments{i,1}(end,6) + segments{i,1}(end,7) - segments{i,1}(1,6);
    x0 = segments{i,1}(1,1);
    x1 = segments{i,1}(end,1) + segments{i,1}(end,2);
    % normalized x points for output curve reading
    segments{i,4} = lin_interpolation(segments{i,1}(:,1) + 0.5*segments{i,1}(:,2), ...
        0, x0, 10, x1);
end

% compute the quadratic coefficients that approximate the contour curve for
% each expertDB segment

quadCoef = zeros(S,3);
for i = 1:S
    x = expertDB{i,1};
    x0 = x(1,1);
    x1 = x(end,1)+x(end,2);
    x = x(:,1)+0.5.*x(:,2);
    quadCoef(i,:) = polyfit(lin_interpolation(x, 0, x0, 10, x1), ...
        expertDB{i,1}(:,5) - expertDB{i,4}, 2);
end

clear H tbk x0 x1

%% Nearest neighbor estimates

if method(idxnn) || method(idxnnq)
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
        nn_calc{i,3} = nn_calc{i,1}/size(segments{i,1},1); % normalized mel. dist.
        matchSeg = expertDB(nn_calc{i,2}(1),:); % nearest neighbor
        nn_calc{i,4} = matchSeg{5}; % segment salience
        nn_calc{i,5} = matchSeg{4} - matchSeg{5}; % match piece mean level
        C = quadCoef(nn_calc{i,2}(1),:); % nearest neighbor quad coefficients
        if method(idxnn)
            nn_calc{i,6} = scaleInterpolate(size(segments{i,1},1), matchSeg{1}(:,5) - matchSeg{4}); % contour
        end
        if method(idxnnq)
            nn_calc{i,7} = C(1).*segments{i,4}.^2 + C(2).*segments{i,4} + C(3); % contour
        end
        nn_calc{i,8} = matchSeg{6}; % relative salience
    end

    clear i j matchSeg x x0 x1

    % Estimate the overall piece dynamics as the duration-weighted average of
    % overall piece dynamics of each segment's nearest neighbor
    overall_nn = [nn_calc{:,5}]*vertcat(segments{:,3})./sum([segments{:,3}]);
    
    % prepare output
    
    if method(idxnnq)
        output_nnq = segments;
        output_nnq(:,3) = nn_calc(:,4); % salience
        output_nnq(:,4) = nn_calc(:,7); % contour
        for i = 1:s
            output_nnq{i,1}(:,5) = overall_nn + output_nnq{i,3} + output_nnq{i,4};
        end
        output_nnq(:,5) = nn_calc(:,8); % relative salience
    end
    if method(idxnn)
        output_nn = segments;
        output_nn(:,3) = nn_calc(:,4); % salience
        output_nn(:,4) = nn_calc(:,6); % contour
        for i = 1:s
            output_nn{i,1}(:,5) = overall_nn + output_nn{i,3} + output_nn{i,4};
        end
        output_nn(:,5) = nn_calc(:,8); % relative salience
    end
    
end

%% kNN estimates

if method(idxknn) || method(idxwknn)
    
    % weight calculations 

    % weights follow an inverse multiquadric kernel
    w = 1./(scores+c);
    % some normalization
    %w = w./(diag(std(w,0,2))*ones(size(w)));

    % gaussian kernel
    %w = exp(- scores.^2 ./ c);
    
    % kNN
    w2 = zeros(s,S);
    w3 = zeros(s,S);
    for i = 1:k
        [v, ind] = max(w,[],2);
        w(sub2ind([s,S],1:s,ind')) = 0;
        w3(sub2ind([s,S],1:s,ind')) = v;
        w2(sub2ind([s,S],1:s,ind')) = 1;
    end
    w = w3;
    clear ind w3

    if method(idxknn)
        % Overall dynamics estimate
        overall_knn = w2*(vertcat(expertDB{:,4})-vertcat(expertDB{:,5}))./sum(w2,2);
        overall_knn = [segments{:,3}]*overall_knn./sum([segments{:,3}]);
        
        wCoef = w2 * quadCoef ./ (sum(w2,2)*ones(1,3));
        
        output_knn = segments;
        
        % Salience Estimate
        output_knn(:,3) = num2cell(w2*vertcat(expertDB{:,5})./sum(w2,2));
        % Contour Estimate
        for i = 1:s
            output_knn{i,4} = wCoef(i,1).*segments{i,4}.^2 + ...
                wCoef(i,2).*segments{i,4} + wCoef(i,3);
            output_knn{i,1}(:,5) = overall_knn + output_knn{i,3} + output_knn{i,4};
        end
        % Relative Salience Estimate (in standard deviations)
        output_knn(:,5) = num2cell(w2*vertcat(expertDB{:,6})./sum(w2,2));
        
    end
    
    if method(idxwknn)
        % Overall dynamics estimate
        overall_wknn = w*(vertcat(expertDB{:,4})-vertcat(expertDB{:,5}))./sum(w,2);
        overall_wknn = [segments{:,3}]*overall_wknn./sum([segments{:,3}]);
        
        wCoef = w * quadCoef ./ (sum(w,2)*ones(1,3));
        
        output_wknn = segments;
        
        % Salience Estimate
        output_wknn(:,3) = num2cell(w*vertcat(expertDB{:,5})./sum(w,2));
        % Contour Estimate
        for i = 1:s
            output_wknn{i,4} = wCoef(i,1).*segments{i,4}.^2 + ...
                wCoef(i,2).*segments{i,4} + wCoef(i,3);
            output_wknn{i,1}(:,5) = overall_wknn + output_wknn{i,3} + output_wknn{i,4};
        end
        % Relative Salience Estimate (in standard deviations)
        output_wknn(:,5) = num2cell(w*vertcat(expertDB{:,6})./sum(w,2));
        
    end
end

if sum(method) > 1
    output = output_wknn;
    outputB = output_knn;
    outputC = output_nnq;
    outputD = output_nn;
elseif method(idxwknn)
    output = output_wknn;
elseif method(idxknn)
    output = output_knn;
elseif method(idxnnq)
    output = output_nnq;
else method(idxnn)
    output = output_nn;
end

end