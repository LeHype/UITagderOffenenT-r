function setLegColorAlpha(leg,options)
    arguments
        % mandatory
            leg
        % optional 
            options.LegColor (1,3) double = 1*[1 1 1];
            options.LegAlpha (1,1) double = 0.6;
    end
    struct2CallerWS(options)
    set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[LegColor';LegAlpha]));
end