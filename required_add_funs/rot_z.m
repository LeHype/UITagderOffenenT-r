function xy_mat_rot = rot_z(xy_mat,alpha,unit)
if nargin<3
    unit='rad';
end
switch unit
    case 'rad'
        a = alpha;
    case 'deg'
        a = alpha/180*pi;
    otherwise
        error('Choose either "rad" or "deg"')
end
rot_mat = [ cos(a) -sin(a)
            sin(a)  cos(a) ];
        
xy_mat_rot = rot_mat*xy_mat;

end