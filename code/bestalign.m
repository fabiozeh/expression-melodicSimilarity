% Computes the best alignment of two numeric sequences according
% to the Needleman-Wunsch algorithm with canonical edit distance
% weights for simple insertion and deletion. The function outputs
% a cell array 'scores' such that:
% scores{1} is the computed length(seq1)+1:length(seq2)+1 matrix
% of scores where end:end is the edit distance between the sequences.
% scores{2} is one possible best alignment path with the precedence
% error < insertion < deletion (of elements in seq2 when compared
% to seq1). scores{2}(a,b,d) is the index of the d'th dimension of
% of the previous cell in scores{1} matrix of the computed path.
% (i.e.: since scores{1}(end:end) is the last cell in the path, 
% scores{1}(scores{2}(end,end,1),scores{2}(end,end,2)) is the
% second-to-last.)
function scores = bestalign(seq1, seq2)

cost_m = 0;
cost_e = 1;
cost_i = 1;
cost_d = 1;

if nargin < 3
    alignseq2 = 0;
end;

path = ones(length(seq1)+1,length(seq2)+1,2);
score = zeros(length(seq1)+1,length(seq2)+1);
score(1,:) = 0:length(seq2);
score(:,1) = 0:length(seq1);
path(1,2:end,2) = 1:length(seq2);
path(2:end,1,1) = 1:length(seq1);
for i = 1:length(seq1)
    for j = 1:length(seq2)
        if (seq1(i)==seq2(j)), match_score = cost_m; else match_score = cost_e; end;
        m_e = score(i,j) + match_score; % match or error
        i_1 = score(i,j+1) + cost_i; % indel 1
        i_2 = score(i+1,j) + cost_d; % indel 2
        if (m_e <= i_1)
            if (m_e <= i_2)
                score(i+1, j+1) = m_e;
                path(i+1, j+1, :) = [i j];
            else
                score(i+1, j+1) = i_2;
                path(i+1, j+1, :) = [i+1 j];
            end
        elseif (i_1 <= i_2)
            score(i+1, j+1) = i_1;
            path(i+1, j+1, :) = [i j+1];
        else
            score(i+1, j+1) = i_2;
            path(i+1, j+1, :) = [i+1 j];
        end
    end
end
scores = {score path};
