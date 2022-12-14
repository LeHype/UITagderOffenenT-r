function h = circle3(pos,d,options)        
    arguments
        % mandatory          
            pos (1,2) double
            d   (1,1) double
        % optional 
            options.FaceColor = 'blue'
            options.EdgeColor = 'black'
            options.FaceAlpha = 1
            options.LineWidth = 1.5
            options.LineStyle = '-'
            options.ScaleX    = 1
            options.ScaleY    = 1
            options.nEdges    = 500
    end
    struct2CallerWS(options)
    
    x = pos(1);
    y = pos(2);
    theta = linspace(0,2*pi,nEdges);
    xunit = ScaleX*d/2*cos(theta) + x;
    yunit = ScaleY*d/2*sin(theta) + y;
    v = [xunit;yunit]';
    h = patch('Faces',1:nEdges,'Vertices',v,...
              'FaceColor',FaceColor,'EdgeColor',EdgeColor,'FaceAlpha',FaceAlpha,'LineWidth',LineWidth,'Linestyle',LineStyle);
end