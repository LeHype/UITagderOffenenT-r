function SPRK_check(a,b,A,B)
s=length(b);
e=[];
for i=1:s
    for j=1:s
        e = [e; b(i)*A(i,j)+B(j)*a(j,i)-b(i)*B(j)];
    end
end
if sum(abs(e))<length(e)*10^-12
    fprintf('The two Butcher tableaus form a SPRK method\n')
else
    warning('The two Butcher tableaus DO NOT form a SPRK method')
end
end