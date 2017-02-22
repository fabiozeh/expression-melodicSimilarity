function yOut = scaleInterpolate(outputSize, yIn)

if ~iscolumn(yIn), yIn = yIn'; end
xIn = (1:size(yIn,1))';

step = (xIn(end)-1)/(outputSize-1);
xOut = 1 + (0:(outputSize-1)).*step;

function y = interpolatePoint(x)
    p1 = floor(x);
    if p1 + 1 <= size(yIn,1)
        y = yIn(p1) * (p1 + 1.0 - x) + yIn(p1+1) * (x - p1);
    else
        y = yIn(p1);
    end
end

yOut = arrayfun(@(arg)interpolatePoint(arg), xOut);

end