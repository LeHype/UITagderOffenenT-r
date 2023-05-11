function R = rot3D(alpha,RotAxis,options)
    arguments
        % mandatory
            alpha
            RotAxis {mustBeMember(RotAxis,{'x','y','z'})}
        % optional 
            options.Unit {mustBeMember(options.Unit,{'deg','rad'})} = 'rad';
            options.LinDeltaVar = []
    end

    Unit = options.Unit;
    LinDeltaVar = options.LinDeltaVar;
    
    if strcmp(Unit,'deg')==1
        alpha = alpha/180*pi;
    end

    if ~isempty(LinDeltaVar)
        d_alpha = LinDeltaVar;
        sin_ = sin(alpha)+cos(alpha)*d_alpha;
        cos_ = cos(alpha)-sin(alpha)*d_alpha;
    else
        sin_ = sin(alpha);
        cos_ = cos(alpha);
    end

    switch RotAxis
        case 'x'
            R = [ 1       0      0
                  0       cos_  -sin_
                  0       sin_   cos_ ];
        case 'y'
            R = [ cos_    0      sin_
                  0       1      0
                  -sin_   0      cos_ ];
        case 'z'
            R = [ cos_   -sin_   0
                  sin_    cos_   0
                  0       0      1    ];
    end
end