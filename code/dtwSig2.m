function [H, tbk] = dtwSig2(nmat1, nmat2, pitchW, durW, onsetW, contW, fig)
% Weights: from 0 to 1

inv = 'no';

%% Create matrix

H=zeros(size(nmat2,1)+1,size(nmat1,1)+1);
H(1,1)=0;
tbk=[];
for i=2:size(nmat2,1)+1
    H(i,1)=inf;
end
for j=2:size(nmat1,1)+1
    H(1,j)=inf;
end

% calculate pitch contours
cont1 = [0; diff(nmat1(:,4))];
cont2 = [0; diff(nmat2(:,4))];

%put nmat1 and nmat2 in the same octave
if (round(mean(nmat1(:,4))-mean(nmat2(:,4)))~=0)
    if mean(nmat1(:,4))<mean(nmat2(:,4))
        nmat2=shift(nmat2,'pitch',-abs(round(mean(nmat1(:,4))-mean(nmat2(:,4)))));
    else
        nmat1=shift(nmat1,'pitch',-abs(round(mean(nmat1(:,4))-mean(nmat2(:,4)))));
    end
end

for j=2:size(nmat1,1)+1
    for i=2:size(nmat2,1)+1
        costCont = contW*(cont2(i-1)-cont1(j-1))^2; %pitch contour difference mult. by weight factor
        costPitch = pitchW*(nmat2(i-1,4)-nmat1(j-1,4))^2; %pitch difference multiplied by a weight factor
        costDur = durW*(nmat2(i-1,2)-nmat1(j-1,2))^2;%duration difference multiplied by a weight factor
        costOnset = onsetW*(nmat2(i-1,1)-nmat1(j-1,1))^2;%onset difference multiplied by a weight factor
        cost = costCont+costPitch+costDur+costOnset;%total cost
        H(i,j) = cost+min([H(i-1,j),H(i,j-1),H(i-1,j-1)]);
    end
end

%% Trace back

i= size(nmat2,1)+1;% go to top left corner
j= size(nmat1,1)+1;% go to top left corner

tbk=[i,j];%initilize trace back
minDir=1;%inizialize minumun path value
while  (minDir~=0 && minDir~=inf)
    minDir=min ([H(i-1,j),H(i-1,j-1),H(i,j-1)]) ;%find minumun path value
    switch minDir % search for minimum value from top, top left and left cells
        case H(i-1,j) %if top cell is less then
          tbk=[tbk;i-1,j];% (insertion)
          i=i-1;%go to previous top cell
        case H(i-1,j-1)%match
          tbk=[tbk;i-1,j-1];
          i=i-1;%go to previous top left cell
          j=j-1;  
        case H(i,j-1)% deletion
             tbk=[tbk;i,j-1];
             j=j-1;%go back to left cell
     end
end

%% Invert correlated notes
  if strcmp(inv,'yes')%if inverted, correlated inversion must be done this way
      tbk(:,1)=(max(tbk(:,1))-tbk(:,1)+1)+1;%max-num+1... ej.: 12233345 -> 54333221(+1 cause indexes start from 2)
      tbk(:,2)=(max(tbk(:,2))-tbk(:,2)+1)+1;%max-num+1... ej.: 12233345 -> 54333221
  else%else just flip up down the vector
      tbk=flipud(tbk);
  end
 
 %% Complete gaps in trace back vector
 % Notes previous to first aligned couple are considered to be part of the
 % first note
 if tbk(1,1)~=2%if trace back of first position is different to 2 (first note) (i dimension)
     for  i=tbk(1,1)-1:-1:2 %from starting value on first sequence
         tbk=[i,tbk(1,2);tbk]; %fill out initial notes with first corresponcence notes couple
     end
 else if tbk(1,2)~=2%if trace back of first position is different to 2 (first note)(j dimension)
        for  i=tbk(1,2)-1:-1:2 %from starting value on first sequence
             tbk=[tbk(1,1),i;tbk]; %fill out initial notes with first corresponcence notes couple
        end
     end
 end
    
 %similar thing to final notes...
  if tbk(end,1)~=size(H,1)%if trace back of first position is different to the lenght of i dimension
     for  i=tbk(end,1)+1:size(H,1) %from starting value on first sequence to end
         tbk=[tbk;i,tbk(end,2)]; %fill out initial notes with first corresponcence notes couple
     end
 else if tbk(end,2)~=size(H,2)%if trace back of first position is different to 2 (first note)(j dimension)
        for  i=tbk(end,2)+1:size(H,2) %from starting value on second sequence to end
             tbk=[tbk;tbk(end,1),i]; %fill out initial notes with first corresponcence notes couple
        end
     end
  end
 
  %gap filling (on first sequence): do the same for second sequence)!
  i=1;
  while i<=length(tbk)-1
      tbkDif= tbk(i+1,1)-tbk(i,1);
      if tbkDif>1
        for j=1:tbkDif-1
            tbk=[tbk(1:i,:);tbk(i,1)+1,tbk(i,2);tbk(i+1:end,:)];
            i=i+1;
        end
      end
      i=i+1;
  end
  %gap filling (on second sequence)
    i=1;
  while i<=length(tbk)-1
      tbkDif= tbk(i+1,2)-tbk(i,2);
      if tbkDif>1
        for j=1:tbkDif-1
            tbk=[tbk(:,1:i);tbk(i,1),tbk(i,2)+1;tbk(:,i+1:end)];
            i=i+1;
        end
      end
      i=i+1;
  end
%  tbk=tbk-1;%reduce one index position as we used a initial row and column of zeros to create H matrix
  %Invert tbk
  
%remove first coumn (first position of H matrix set to zero) and reduce 1 index as H is zero indexed
tbk=tbk(2:end,:)-1;

if strcmp(fig,'yes')
 %% Plot aligment
 fig1=figure (3);
 image(sqrt(H(2:end,2:end)), 'CDataMapping','scaled');
 set(gca,'YDir','normal');
 set(gca,'FontSize',12);
 figure(gcf);
 hold on
 plot (tbk(:,2),tbk(:,1),'k*')
 xlabel('Score note number','FontSize',14);
 ylabel('Performed note number','FontSize',14);
 hndl=colorbar;
 ylabel(hndl,'Cost','FontSize',14);
 
 %plot piano roll
 pnrll=aligmentPlot(nmat1,nmat2,tbk,1);
end

end

function pnrll=aligmentPlot(nmat,nmat2,tbk,octShift)
%This function plots a piano roll of two songs one octave appart, and draw
%lines betwen correspondig notes, base on the aligment vector tbk. Input
%variables are:
%nmat: midi matrix of score (midi format based on midi toolbox [ref])
%namt: midi matrix of performance
%tbk: Aligment between notes of score and performance
%octShift: octave shift to plot betwen performance and score

nmat=shift(nmat,'pitch',octShift*12);
nmat2=shift(nmat2,'pitch',-octShift*12);
%subplot(2,1,1)
scrsz = get(0,'ScreenSize');% get screen size
%figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2])

pnrll=figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)/2]);%plot piano roll in half of the screen

pianoroll(nmat);
%subplot(2,1,2)
pianoroll(nmat2, 'g', 'hold','num','beat');
hold on;
 %find pairs
%         onset(beat)         + half duration in beats (so marking will be at the middle of note box)
all_x=[(nmat2(tbk(:,1),1)'+nmat2(tbk(:,1),2)'/2); ...
           (nmat(tbk(:,2),1)' +nmat(tbk(:,2),2)'/2)];        %first set of ponts (notes) from second matrix (performed)
       
all_y=[nmat2(tbk(:,1),4)'; ... 
           nmat(tbk(:,2),4)'];%second set of ponts (notes) from second matrix (score)


% % % x = [0 1 1 0; ...
% % %      1 1 0 0];
% % % y = [0 0 1 1; ...
% % %      0 1 1 0];
% % % plot(x,y);
% % % This will plot each line in a different color. To plot all of the lines as black, do this:
% % % 
plot(all_x,all_y);

%% dysplay note numbers
for i=1:size(nmat2,1)
    text(nmat2(i,1),nmat2(i,4)+1,num2str(i));
end

for i=1:size(nmat,1)
    text(nmat(i,1),nmat(i,4)+1,num2str(i));
end

end
