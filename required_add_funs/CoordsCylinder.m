function coords = CoordsCylinder(r,H,N,options)
    arguments
        % mandatory
            r
            H
            N
        % optional 
            options.Translation = zeros(3,1);
            options.Rotation    = eye(3);
            options.RotCenterShift = zeros(3,1);
    end
    struct2CallerWS(options)


    theta_grid = linspace(0,2*pi,N);
    x = r*cos(theta_grid);
    y = r*sin(theta_grid);

    h = H/2;

    X = [x,x];
    Y = [y,y];
    Z = [h*ones(1,N),-h*ones(1,N)];
    
    coords = [X;Y;Z];

    coords = Translation+(Rotation*(coords-RotCenterShift));
end