% Velocities calculated from dBFS values according to dannenberg's 
% "square-law" that says the square-root of RMS values measured from 
% synthesized sounds are linearly correlated to input velocity.
% The lo and hi variables set dB levels for velocities 1 and 120;
function midiout = dbfs2vel_sqrt(midiin)

lo = -36;
hi = -6;

x0 = 10^(lo/40);
x1 = 10^(hi/40);

if size(midiin, 2) > 1

    midiout = midiin;
    
    midiout(:,5) = min(horzcat(127 + zeros(size(midiin,1),1), ...
      max(horzcat(10 + zeros(size(midiin,1),1), ...
      lin_interpolation(10.^(midiin(:,5)./40), 10, x0, 120, x1)), [], 2)), [], 2);

else
    midiout = min(horzcat(127 + zeros(size(midiin,1),1), ...
      max(horzcat(10 + zeros(size(midiin,1),1), ...
      lin_interpolation(10.^(midiin(:,1)./40), 10, x0, 120, x1)), [], 2)), [], 2);
    
end