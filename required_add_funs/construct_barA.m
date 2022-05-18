function barA = construct_barA(A,b)
s = length(b);
barA = NaN*ones(s,s);
for ii=1:s
    for jj=1:s
        barA(ii,jj) = b(jj)*(1-A(jj,ii)/b(ii));
    end
end
end