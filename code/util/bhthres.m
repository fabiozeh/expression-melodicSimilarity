function t = bhthres(v, init)
sum1 = 0;
n1 = 0;
sum2 = 0;
n2 = 0;
if nargin < 2, init = mean(v); end
t = init;
while(1)
    for ii = 1:numel(v)
        if v(ii) < t
            sum1 = sum1 + v(ii);
            n1 = n1 + 1;
        else
            sum2 = sum2 + v(ii);
            n2 = n2 + 1;
        end
    end
    new = (sum1/n1 + sum2/n2)/2;
    if new == t, break; end
    t = new;
    sum1 = 0;
    sum2 = 0;
    n1 = 0;
    n2 = 0;
end
