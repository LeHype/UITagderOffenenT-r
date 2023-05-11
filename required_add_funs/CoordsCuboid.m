function coords = CoordsCuboid(L,D,H,options)
    arguments
        % mandatory
            L
            D
            H
        % optional 
            options.Translation = zeros(3,1);
            options.Rotation    = eye(3);
            options.RotCenterShift = zeros(3,1);
    end
    struct2CallerWS(options)


    l = L/2;
    d = D/2;
    h = H/2;
    
    X = l*[-1  1  1 -1 -1  1  1 -1];
    Y = d*[-1 -1  1  1 -1 -1  1  1];
    Z = h*[-1 -1 -1 -1  1  1  1  1];
    coords = [X;Y;Z];

    coords = Translation+(Rotation*(coords-RotCenterShift));
end