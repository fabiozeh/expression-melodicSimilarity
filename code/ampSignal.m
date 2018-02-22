function asig = ampSignal(nmat, sRate) 

% if necessary, convert velocity values to dbfs
if nmat(1,5) > 0, nmat = vel2dbfs(nmat); end

asig(1:ceil((nmat(end, 6) + nmat(end, 7)) * sRate),1) = -150;
for i = 1:size(nmat,1)

    if (nmat(i,4) == 0), a = -120; else a = nmat(i,5); end
    st = 1+round(nmat(i,6)*sRate);
    e = st+round(nmat(i,7)*sRate);
    asig(st:e) = a;
    
end

end