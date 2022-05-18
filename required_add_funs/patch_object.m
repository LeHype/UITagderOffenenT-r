function all_patches = patch_object(obj,obj_colors,FaceAlpha)
    for ii=1:length(obj)
        part_ii = obj{ii};
        all_patches(ii) = patch(part_ii(1,:),part_ii(2,:),obj_colors{ii},'FaceAlpha',FaceAlpha);
    end
end