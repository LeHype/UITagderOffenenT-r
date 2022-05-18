function inds_keep = closestPointsInCurve(curve_xy)    
    n_scan    = ceil(size(curve_xy,2)/3.5);
    n_neighb  = ceil(n_scan/3);
    max_min_l = norm(curve_xy(:,1)-curve_xy(:,2))*10;    

    info_mat = [];
    for jj=1:length(curve_xy)
        ints = [max([1,(jj-n_scan)]):(jj-1-n_neighb), (jj+1+n_neighb):min([length(curve_xy),(jj+n_scan)])];
        dist_mat = curve_xy(:,ints) - curve_xy(:,jj);
        l_vec = vecnorm(dist_mat,2,1);
        [min_l,ind_min_l] = min(l_vec);
        if min_l<max_min_l
            info_mat = [info_mat, [min_l;jj;ints(ind_min_l)]];
        end
    end
%     info_mat(2,:)
    
    try
        inds_help = find(diff(info_mat(2,:))~=1);
        n_ints = length(inds_help)+1;
        ind_se_info_mat = zeros(2,n_ints);
        for ii=1:n_ints
            if ii==1
                is = 1;
            else
                is = inds_help(ii-1)+1;
            end
            if ii==n_ints
                ie = size(info_mat,2);
            else
                ie = inds_help(ii);
            end
            ind_se_info_mat(:,ii) = [is;ie];
        end
        
        info_mat_red = [];
        for ii=1:n_ints
            info_mat_slice = info_mat(:,ind_se_info_mat(1,ii):ind_se_info_mat(2,ii));
            [~,ind_min_l] = min(info_mat_slice(1,:));
            info_mat_red = [info_mat_red, info_mat_slice(:,ind_min_l)];
            if info_mat_red(3,end)<info_mat_red(2,end)
                info_mat_red([2,3],end) = info_mat_red([3,2],end);
            end
        end
            
        ind_delete = [];
        for ii=1:size(info_mat_red,2)
            ind_delete = [ind_delete, info_mat_red(2,ii):info_mat_red(3,ii)];
        end
    
        inds_keep = setdiff(1:size(curve_xy,2),ind_delete);
    catch
        inds_keep = 1:size(curve_xy,2);
    end


end
