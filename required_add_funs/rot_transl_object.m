function obj_out = rot_transl_object(obj_in,angle,angle_unit,transl,rot_axis,indiv_part)
    if isempty(rot_axis)
        rot_axis = [0 0]';
    end
    obj_out = obj_in;
    for ii=1:length(obj_in)
        part_ii = obj_in{ii};
        if ~isempty(indiv_part) && any(ii==indiv_part.nr)
            part_ii = rot_transl_part(part_ii,indiv_part.angle(ii==indiv_part.nr),angle_unit,[0;0],indiv_part.axis(:,ii==indiv_part.nr));
        end
        part_ii = rot_transl_part(part_ii,angle,angle_unit,transl,rot_axis);

        obj_out(ii) = {part_ii};
    end
end