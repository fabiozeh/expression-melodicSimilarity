function dumpMidis(folder)

if isunix(), sep = '/'; else sep = '\'; end

d = dir(folder);
d = d(contains({d.name},'.wav'));
pieceList = {d.name};
pieceList = strcat([folder sep], pieceList);
for piece = pieceList
    [feats, fData] = createExpertDB(piece, 0, 0);
    performance = vertcat(feats{:,1});
    perfmid = performance(:,1:7);
    perfmid(:,6:7) = performance(:,9:10);
    deadpmid = performance(:,1:7);
    deadpmid(:,5) = fData(1,1);
    name = piece{1};
    writemidi(deadpmid, [name(1:end-3) '_deadpan.mid']);
    writemidi(perfmid, [name(1:end-3) '_performance.mid']);
end