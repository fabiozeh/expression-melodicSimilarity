%% Outputs a vector of size outputSize by linearly interpolating values from yIn
function yOut = scaleInterpolate(outputSize, yIn)

if ~iscolumn(yIn), yIn = yIn'; end
xIn = (1:size(yIn,1))';

step = (xIn(end)-xIn(1))/(outputSize-1);
xOut = (1 + (0:(outputSize-1)).*step)';
yIn = [yIn; yIn(end)];
p1 = floor(xOut);
yOut = yIn(p1) .* (p1 + 1.0 - xOut) + yIn(p1+1) .* (xOut - p1);

end