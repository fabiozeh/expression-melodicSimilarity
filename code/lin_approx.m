function y = lin_approx(newx, x, fx)

y = zeros(size(newx,1),1);
for i = 1:size(newx,1)
x0 = max(x(x <= newx(i)));
if isempty(x0), x0 = x(1); end
y0 = fx(x == x0);
if size(y0,1) > 1, y0 = mean(y0); end
x1 = min(x(x >= newx(i)));
if isempty(x1), x1 = x(end); end
y1 = fx(x == x1);
if size(y1,1) > 1, y1 = mean(y1); end

if x0 == x1
    y(i) = y0;
else
    b = (y1 - y0)./(x1 - x0);
    a = y0 - b.*x0;
    y(i) = a + b.*newx(i);
end
end

end