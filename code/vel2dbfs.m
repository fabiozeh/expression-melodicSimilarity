function midiout = vel2dbfs(midiin)

x0 = 80;
y0 = -20; % mean dbfs
x1 = 24; %x0 - 60; % minimum audible dbfs
y1 = -70; % 20% midi velocity

midiout = midiin;
midiout(:,5) = lin_interpolation(midiout(:,5), y0, x0, y1, x1);
midiout(midiout(:,5) < -70,5) = lin_interpolation(midiout(midiout(:,5) < -70,5), -80, min(midiout(:,5)), -70, -70);
midiout(midiout(:,5) > -15,5) = lin_interpolation(midiout(midiout(:,5) > -15,5), -4, max(midiout(:,5)), -15, -15);

end