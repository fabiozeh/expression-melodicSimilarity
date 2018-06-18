% Velocities calculated from dBFS values according to dannenberg's 
% "square-law" that says the square-root of RMS values measured from 
% synthesized sounds are linearly correlated to input velocity.
% The lo and hi variables set dB levels for velocities 1 and 120;
function midiout = dbfs2vel_sqrt(midiin)

midiout = midiin;

lo = -38;
hi = -6;

x0 = 10^(lo/40);
x1 = 10^(hi/40);
midiout(:,5) = lin_interpolation(10.^(midiin(:,5)./40), 1, x0, 120, x1);

end