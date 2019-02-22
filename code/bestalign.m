% Computes the best alignment of two numeric sequences according
% to the Needleman-Wunsch algorithm with canonical edit distance
% weights for simple insertion and deletion. The function outputs:
% score, which is the computed length(seq1)+1:length(seq2)+1 matrix
% of scores where end:end is the edit distance between the sequences;
% and path, which is one possible best alignment path with the precedence
% error < insertion < deletion (of elements in seq2 when compared to seq1).
function [score, path] = bestalign(seq1, seq2)

cost_m = 0;
cost_e = 1;
cost_i = 1;
cost_d = 1;

score = zeros(length(seq1)+1,length(seq2)+1);
score(1,:) = 0:length(seq2);
score(:,1) = 0:length(seq1);
for i = 1:length(seq1)
    for j = 1:length(seq2)
        if (seq1(i)==seq2(j)), match_score = cost_m; else match_score = cost_e; end;
        m_e = score(i,j) + match_score; % match or error
        i_1 = score(i,j+1) + cost_i; % indel 1
        i_2 = score(i+1,j) + cost_d; % indel 2
        if (m_e <= i_1)
            if (m_e <= i_2)
                score(i+1, j+1) = m_e;
            else
                score(i+1, j+1) = i_2;
            end
        elseif (i_1 <= i_2)
            score(i+1, j+1) = i_1;
        else
            score(i+1, j+1) = i_2;
        end
    end
end
k = 2;
path(1,:) = [i, j];
while (i > 1 && j > 1)
    d = score(i-1,j-1);
    u = score(i,j-1);
    l = score(i-1,j);
    if (d <= u)
        if (d <= l)
           path(k,:) = [i-1; j-1];
        else
           path(k,:) = [i-1; j];
        end
    else
        if (u <= l)
            path(k,:) = [i; j-1];
        else
            path(k,:) = [i-1; j];
        end
    end
    i = path(k,1);
    j = path(k,2);
    k = k+1;
end
path = path(end:-1:1,:);
end
