sf1 = cell2mat(sound_file_1);
sf2 = cell2mat(sound_file_2);

%cond = sound_sample > 0; % first chart
cond = ismember(sound_sample, [31, 61, 71]); % second chart (filtered)

acMask = sf1(cond) == 'a' & sf2(cond) == 'c';
caMask = sf1(cond) == 'c' & sf2(cond) == 'a';
abMask = sf1(cond) == 'a' & sf2(cond) == 'b';
baMask = sf1(cond) == 'b' & sf2(cond) == 'a';
bcMask = sf1(cond) == 'b' & sf2(cond) == 'c';
cbMask = sf1(cond) == 'c' & sf2(cond) == 'b';

h = human(cond);
p = preferred(cond);

aoch = sum((acMask & h == 1) | (caMask & h == 2));
aocp = sum((acMask & p == 1) | (caMask & p == 2));
ac = sum(acMask | caMask);
aobh = sum((abMask & h == 1) | (baMask & h == 2));
aobp = sum((abMask & p == 1) | (baMask & p == 2));
ab = sum(abMask | baMask);
boch = sum((bcMask & h == 1) | (cbMask & h == 2));
bocp = sum((bcMask & p == 1) | (cbMask & p == 2));
bc = sum(bcMask | cbMask);

figure
barfig = bar([aoch aocp ac; aobh aobp ab; boch bocp bc]);
%title('Aggregate results of perceptual survey'); % first graph
title('Filtered results of perceptual survey'); % second graph
set(gca,'xticklabel',{'Performer over Baseline','Performer over Model',...
    'Model over Baseline'});
legend('Human-like', 'Preferred', 'Total');



% fm = ismember(user_id,id1(if_lessons_years>0|hours_practice>0));
% meandis = zeros(24,2);
% control = zeros(24,3);
% ind = 1;
% for i = 1:8
%     for j = 1:3
%         smpl = 10*i+j;
%         meandis(ind,:) = [smpl,mean(distinction(sound_sample(fm) == smpl))];
%         control(ind,1:2) = [smpl,sum(((sf1(fm)=='a' & sf2(fm)=='c') | ...
%             (sf1(fm)=='c' & sf2(fm)=='a')) & sound_sample(fm) == smpl)];
%         control(ind,3) = sum(((sf1(fm)=='a' & sf2(fm)=='c' & human(fm) == 1) | ...
%             (sf1(fm)=='c' & sf2(fm)=='a' & human(fm) ==2)) & sound_sample(fm) == smpl).*1.0/ ...
%             control(ind,2);
%         ind = ind+1;
%     end
% end
% 
% sum(((sf1(fm)=='c' & sf2(fm)=='a' & human(fm) == 2) | (sf2(fm)=='c' & sf1(fm)=='a' & human(fm) == 1)) & ismember(sound_sample(fm),meandis(meandis(:,2)>=3,1)))