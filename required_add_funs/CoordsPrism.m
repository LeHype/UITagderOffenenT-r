function coords = CoordsPrism(G,D,H,dx,options)
    arguments
        % mandatory
            G
            D
            H
        % optional 
            dx = 0;
            options.Translation = zeros(3,1);
            options.Rotation    = eye(3);
            options.RotCenterShift = zeros(3,1);
    end
    struct2CallerWS(options)

    g = G/2;
    d = D/2;
    
    X = -g*[-1  1  0  -1  1  0] ...
         + [0   0  dx  0  0  dx];
    Y = d*[-1 -1 -1 1 1 1];
    Z = H*[0 0 1 0 0 1];
    coords = [X;Y;Z];

    coords = Translation+(Rotation*(coords-RotCenterShift));
end