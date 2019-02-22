% Perform dynamics estimation for a given score.
% Parameters:
%  score: nmat of notes in the style of miditoolbox or cell array of
% phrases as output by findPhrases function.
%  L: the desired mean velocity level (dbFS)
%  R: the desired dynamic range (dbFS)
%  expertDB: the training set, created with createExpertDB.m
%  method: 'w-knn', 'knn', 'nn-quad', 'nn', or 'all'
%  k: number of neighbors for k-nn or w-knn (typical: 2-10). Can be an array
% for multiple estimations.
%  c: RBF kernel parameter for w-knn (typically < 1e-3 but should be
%    cross-validated for optimal results.
% Output:
%  if method = 'all', an array of 4 cell arrays with output data, otherwise
% a single cell array with the desired method. If method includes knn or 
% w-knn, there will be as many 3rd dimensions in the corresponding cell
% arrays as k values in input. Content of the cell array is as follows.
%  rows: 1 row per phrase, as in input or detected with findPhrases.m
%  column 1: score nmat with modeled velocity levels.
%  column 2: index of first phrase note in piece.
%  column 3: alpha prediction (dynamic salience)
%  column 4: beta prediction (dynamic range)
%  column 5: gamma prediction (dynamic contour)
%  column 6: parametric gama prediction (parabolic contour coefficients)
%  column 7: mean melodic distance to chosen training set instances
function [output, outputB, outputC, outputD] = dynamicsEstimation(score, L, R, expertDB, method, k, c)

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

if (method(idxknn) || method(idxwknn)) && nargin < 6, k = 7; c = 1e-9; end
if (method(idxknn) || method(idxwknn)) && nargin < 7, c = 1e-9; end

%% segment the score into melodic phrases, if needed
if ~iscell(score)
    segments = findPhrases(score);
else
    segments = score(:,1:2);
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
expertDBxs = cell(S,1);
for i = 1:S
    x = expertDB{i,1};
    x0 = x(1,1);
    x1 = x(end,1)+x(end,2);
    x = x(:,1)+0.5.*x(:,2);
    expertDBxs{i} = lin_interpolation(x, 0, x0, 10, x1);
    quadCoef(i,:) = polyfit(expertDBxs{i}, ...
        expertDB{i,10},2);
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
        nn_calc{i,4} = matchSeg{8}; % alpha (segment salience coeff.)
        nn_calc{i,5} = matchSeg{9}; % beta (segment std. coeff.)
        C = quadCoef(nn_calc{i,2}(1),:); % nearest neighbor quad coefficients
        nn_calc{i,8} = C;
        if method(idxnn)
            nn_calc{i,6} = lin_approx(segments{i,4}, ...
                expertDBxs{nn_calc{i,2}(1)}, matchSeg{10});% contour
        end
        if method(idxnnq)
            nn_calc{i,7} = C(1).*segments{i,4}.^2 + C(2).*segments{i,4} + C(3); % contour
        end
        
    end

    clear i j matchSeg x x0 x1

    % prepare output
    
    if method(idxnnq)
        output_nnq = segments;
        output_nnq(:,3) = nn_calc(:,4); % alpha (salience coeff.)
        output_nnq(:,4) = nn_calc(:,5); % beta (std. coeff.)
        output_nnq(:,5) = nn_calc(:,7); % gamma (note contour coeff.)
        output_nnq(:,6) = nn_calc(:,8); % parametric gamma (parabolic) coeff.
        output_nnq(:,7) = nn_calc(:,3); % normalized melodic distance
        for i = 1:s
            output_nnq{i,1}(:,5) = L + R.*(output_nnq{i,3} + ...
                output_nnq{i,4}.*output_nnq{i,5});
            output_nnq(i,1) = {dbfs2vel_sqrt(output_nnq{i,1})};
        end
        
    end
    if method(idxnn)
        output_nn = segments;
        output_nn(:,3) = nn_calc(:,4); % alpha (salience coeff.)
        output_nn(:,4) = nn_calc(:,5); % beta (std. coeff.)
        output_nn(:,5) = nn_calc(:,6); % gamma (note contour coeff.)
        %output_nn(:,6) = []; % no parametric gamma for you! 
        output_nn(:,7) = nn_calc(:,3); % normalized melodic distance
        for i = 1:s
            output_nn{i,1}(:,5) = L + R.*(output_nn{i,3} + ...
                output_nn{i,4}.*output_nn{i,5});
            output_nn(i,1) = {dbfs2vel_sqrt(output_nn{i,1})};
        end
        
    end
    
end

%% kNN estimates

if method(idxknn) || method(idxwknn)
    
    output_knn = cell(s,7,length(k));
    output_wknn = cell(s,7,length(k));
    for kind = 1:length(k)
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
        avgDist = zeros(s,1);
        for i = 1:(k(kind))
            [v, ind] = max(w,[],2);
            w(sub2ind([s,S],1:s,ind')) = 0;
            w3(sub2ind([s,S],1:s,ind')) = v;
            w2(sub2ind([s,S],1:s,ind')) = 1;
            avgDist = avgDist + v;
        end
        w = w3;
        avgDist = avgDist ./ k(kind);
        
        clear ind w3

        if method(idxknn)

            wCoef = w2 * quadCoef ./ (sum(w2,2)*ones(1,3));

            output_knn(:,1:2,kind) = segments(:,1:2);

            % Alpha Estimate
            output_knn(:,3,kind) = num2cell(w2*vertcat(expertDB{:,8})./sum(w2,2));
            % Beta estimate
            output_knn(:,4,kind) = num2cell(w2*vertcat(expertDB{:,9})./sum(w2,2));
            % Gamma Estimate
            for i = 1:s
                output_knn{i,5,kind} = wCoef(i,1).*segments{i,4}.^2 + ...
                    wCoef(i,2).*segments{i,4} + wCoef(i,3);
                output_knn{i,1,kind}(:,5) = L + R.*(output_knn{i,3,kind} + ...
                    output_knn{i,4,kind}.*output_knn{i,5,kind});
                output_knn(i,1,kind) = {dbfs2vel_sqrt(output_knn{i,1,kind})};
                output_knn(i,6,kind) = {wCoef(i,:)}; % parametric gamma coefficients
                output_knn(i,7,kind) = {avgDist(i,:)};
            end
            
        end

        if method(idxwknn)

            wCoef = w * quadCoef ./ (sum(w,2)*ones(1,3));

            output_wknn(:,1:2,kind) = segments(:,1:2);

            % Alpha Estimate
            output_wknn(:,3,kind) = num2cell(w*vertcat(expertDB{:,8})./sum(w,2));
            % Beta Estimate
            output_wknn(:,4,kind) = num2cell(w*vertcat(expertDB{:,9})./sum(w,2));
            % Gamma Estimate
            for i = 1:s
                output_wknn{i,5,kind} = wCoef(i,1).*segments{i,4}.^2 + ...
                    wCoef(i,2).*segments{i,4} + wCoef(i,3);
                output_wknn{i,1,kind}(:,5) = L + R.*(output_wknn{i,3,kind} + ...
                    output_wknn{i,4,kind}.*output_wknn{i,5,kind});
                output_wknn(i,1,kind) = {dbfs2vel_sqrt(output_wknn{i,1,kind})};
                output_wknn(i,6,kind) = {wCoef(i,:)}; % parametric gamma coefficients
                output_wknn(i,7,kind) = {avgDist(i,:)};
            end

        end
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
else %method(idxnn)
    output = output_nn;
end

end