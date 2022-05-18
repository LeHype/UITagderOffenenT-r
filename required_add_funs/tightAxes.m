function tightAxes(options)
% based on https://de.mathworks.com/matlabcentral/answers/352024-programmatically-performing-expand-axes-to-fill-figure
    arguments
        % mandatory
            % nothing
        % optional 
            options.AxesH = gca;
            options.PercW = 95/100;
            options.PercH = 97/100;
    end
    struct2CallerWS(options)
    InSet = get(AxesH, 'TightInset');
    
    % before scaling
    l = InSet(1);
    b = InSet(2);
    w = 1-InSet(1)-InSet(3);
    h = 1-InSet(2)-InSet(4);
    
    % scaling
    w = PercW*w;
    l = l+(1-PercW)*w/2;
    h = PercH*h;
    b = b+(1-PercH)*h/2;
    
    set(AxesH, 'Position', [l b w h])
end

