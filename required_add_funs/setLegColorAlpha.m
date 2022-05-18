function setLegColorAlpha(leg,LegColor,LegAlpha)
    arguments
        % mandatory
            leg
        % optional 
            LegColor (1,3) double = 0.95*[1 1 1];
            LegAlpha (1,1) double = 0.6;
    end
    set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[LegColor';LegAlpha]));
end