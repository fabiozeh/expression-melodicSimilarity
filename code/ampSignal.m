function asig = ampSignal(nmat, sRate) 

nmat = vel2dbfs(nmat);
asig = zeros(ceil((nmat(end, 6) + nmat(end, 7)) * sRate), 1);
for i = 1:size(nmat,1)

    if (nmat(i,4) == 0), a = -120; else a = nmat(i,5); end
    st = 1+round(nmat(i,6)*sRate);
    e = st+round(nmat(i,7)*sRate);
    asig(st:e) = a;
    
end

end