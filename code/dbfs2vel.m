function midiout = dbfs2vel(midiin)

minvel = min(midiin(:,5));
maxvel = max(midiin(:,5));
x0 = minvel + (maxvel - minvel) / 2;
y0 = 64;
x1 = x0 - 48;
y1 = 1;

midiout = midiin;
midiout(:,5) = lin_interpolation(midiout(:,5), y0, x0, y1, x1);
midiout(midiout(:,5) < 1,5) = 1;
midiout(midiout(:,5) > 127,5) = 127;

end