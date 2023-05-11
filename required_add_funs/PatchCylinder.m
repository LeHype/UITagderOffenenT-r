function varargout = PatchCylinder(coords,options)
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

    N = size(coords,2)/2;

    faces_help = [1 1+N 2+N 2];
    faces_walls = [faces_help+(0:N-2)';
                   N N+N 1+N 1];
    faces_topbot = [1:N;
                    N+1:2*N];
    faces = [faces_walls,NaN*ones(N,N-4); %https://de.mathworks.com/matlabcentral/answers/176339-how-to-draw-multiple-patches-with-different-number-of-vertices
             faces_topbot];

    p_handle = patch(Axes,'Vertices',vertices,'Faces',faces, ...
                     'FaceColor',FaceColor,'EdgeColor',EdgeColor, ...
                     'LineWidth',LineWidth,'FaceAlpha',FaceAlpha);


    varargout{1} = p_handle;
    varargout{2} = coords;
    varargout = varargout(1:nargout);

end