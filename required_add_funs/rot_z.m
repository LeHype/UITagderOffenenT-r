function xy_mat_rot = rot_z(xy_mat,alpha,unit)

arguments
    %mandatory
        xy_mat
        alpha        
    %optional
        unit {mustBeMember(unit,{'rad','deg'})} = 'rad'
end

switch unit
    case 'rad'
        a = alpha;
    case 'deg'
        a = alpha/180*pi;
end

rot_mat = [ cos(a) -sin(a)
            sin(a)  cos(a) ];
        
xy_mat_rot = rot_mat*xy_mat;

end