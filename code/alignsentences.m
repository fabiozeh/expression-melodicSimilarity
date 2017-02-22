function alignment = alignsentences(sent1,sent2)

M = eye(256).*2 - 1;
[~, path(:,2), path(:,1)] = bioinfoprivate.affinegapmex(uint8(sent1), uint8(sent2), -2, -1, M, 1);
path = path(sum(path,2)>0,:);
path = flipud(path);
alignment = repmat(('--')',1,size(path,1));
alignment(1,path(:,1)>0) = sent2;
alignment(2,path(:,2)>0) = sent1;