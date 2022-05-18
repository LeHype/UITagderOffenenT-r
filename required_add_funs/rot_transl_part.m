function part_out = rot_transl_part(part_in,angle,angle_unit,transl,rot_axis)
    if nargin<5
        rot_axis = [0 0]';
    end
    part_out = transl+(+rot_axis+rot_z(part_in-rot_axis,angle,angle_unit));
end