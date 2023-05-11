function varargout = PatchPrism(coords,options)
    arguments
        % mandatory
            coords
        % optional 
            options.Translation = zeros(3,1);
            options.Rotation    = eye(3);
            options.RotCenterShift = zeros(3,1);
            options.FaceColor   = 'red';
            options.EdgeColor   = 'black';
            options.FaceAlpha   = 0.5;
            options.LineWidth   = 1.2;
            options.Axes        = gca;
    end
    struct2CallerWS(options)

    coords = Translation+(Rotation*(coords-RotCenterShift));

    vertices = coords';
    faces = [1 2 3 NaN;   
             4 5 6 NaN; %	       6          
             2 5 6 3;   %	      /\   
             1 4 6 3];  %	     /  \               
                        %	    4----5
                        %	   /	/
                        %	  /	   /
                        %    3    /       
                        %   /\   /         
                        %  /  \ /          
                        % 1----2

    p_handle = patch(Axes,'Vertices',vertices,'Faces',faces, ...
                     'FaceColor',FaceColor,'EdgeColor',EdgeColor, ...
                     'LineWidth',LineWidth,'FaceAlpha',FaceAlpha);

    varargout{1} = p_handle;
    varargout{2} = coords;
    varargout = varargout(1:nargout);

end