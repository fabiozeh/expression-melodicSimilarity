function midiout = dbfs2vel(midiin)

%minvel = quantile(midiin(:,5), 0.2); %min(midiin(:,5));
%maxvel = quantile(midiin(:,5), 0.8); %max(midiin(:,5));
%x0 = minvel + (maxvel - minvel) / 2; % mean observed dbfs
%y0 = 60; % mean midi velocity
%x1 = minvel; %x0 - 60; % minimum audible dbfs
%y1 = 24; % 20% midi velocity

x0 = -36;
y0 = 1;
x1 = -6;
y1 = 127;

midiout = midiin;
midiout(:,5) = lin_interpolation(midiout(:,5), y0, x0, y1, x1);
%midiout(midiout(:,5) < 24,5) = lin_interpolation(midiout(midiout(:,5) < 24,5), 8, min(midiout(:,5)), 24, 24);
%midiout(midiout(:,5) > 111,5) = lin_interpolation(midiout(midiout(:,5) > 111,5), 127, max(midiout(:,5)), 111, 111);
midiout(midiout(:,5) < 1,5) = deal(1);
midiout(midiout(:,5) > 127,5) = deal(127);
end