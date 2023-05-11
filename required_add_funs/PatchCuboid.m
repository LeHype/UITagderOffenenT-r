function varargout = PatchCuboid(coords,options)
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
    faces = [1 2 6 5;  	%      8----7	
		     2 3 7 6;	%     /|   /|
		     3 4 8 7;	%    / |  / |
		     4 1 5 8;	%   /  5-/--6
		     1 2 3 4;	%  /  / /  /
		     5 6 7 8];	% 4--/-3  /
					    % | /  | /
					    % |/   |/
					    % 1----2
    p_handle = patch(Axes,'Vertices',vertices,'Faces',faces, ...
                     'FaceColor',FaceColor,'EdgeColor',EdgeColor, ...
                     'LineWidth',LineWidth,'FaceAlpha',FaceAlpha);

    varargout{1} = p_handle;
    varargout{2} = coords;
    varargout = varargout(1:nargout);

end