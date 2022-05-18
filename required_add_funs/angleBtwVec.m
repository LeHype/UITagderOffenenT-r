function alpha = angleBtwVec(u,v,CasadiAtan2)
if nargin==3 && CasadiAtan2==1
    u = u([2,1]);
    v = v([2,1]);
end

% x=v(1)=1, y=v(2)=0
alpha = atan2(v(2), v(1)) - atan2(u(2), u(1));
% alpha = atan2Appr(v(2), v(1)) - atan2Appr(u(2), u(1));
% alpha = atan(v(2)/v(1)) - atan(u(2)/u(1));

% alpha = acos((u'*v)/sqrt((u'*u)*(v'*v)));

end