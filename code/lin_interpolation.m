function y = lin_interpolation(x, y0, x0, y1, x1)
    b = (y1 - y0)./(x1 - x0);
    a = y0 - b.*x0;
    y = a + b.*x;
end