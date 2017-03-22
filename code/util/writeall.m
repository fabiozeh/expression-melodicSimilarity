function writeall(segments, folder)
if isunix(), sep = '/'; else sep = '\'; end
for ii = 1:length(segments)
    writemidi(segments{ii,1}(:,1:7),[folder sep 'segm' num2str(ii) '.mid']);
end
clear ii sep
end