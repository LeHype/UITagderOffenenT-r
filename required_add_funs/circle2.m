function h = circle2(pos,d,color)        
    x = pos(1);
    y = pos(2);
    th = linspace(0,2*pi,500);
    xunit = d/2 * cos(th) + x;
    yunit = d/2 * sin(th) + y;
    if nargin<3
        h = plot(xunit, yunit);
    else
        h = plot(xunit, yunit,'Color',color);
    end
end