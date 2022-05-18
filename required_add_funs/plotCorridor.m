function plotCorridor(xy_corridor,max_abs_theta_2)

%% create corridor
h = 5*10^-3;
% theta_1_grid = -1:h:0;
theta_1_grid = 0:h:1;

try
    xy_eval_upper = xy_corridor(theta_1_grid, max_abs_theta_2);
    xy_eval_lower = xy_corridor(theta_1_grid,-max_abs_theta_2);
catch
    xy_eval_upper = repmat(xy_corridor([theta_1_grid(1);  max_abs_theta_2]),1,length(theta_1_grid));
    xy_eval_lower = repmat(xy_corridor([theta_1_grid(1); -max_abs_theta_2]),1,length(theta_1_grid));
    for ii=1:length(theta_1_grid)
        xy_eval_upper(:,ii) = xy_corridor([theta_1_grid(ii);  max_abs_theta_2]);
        xy_eval_lower(:,ii) = xy_corridor([theta_1_grid(ii); -max_abs_theta_2]);
    end
    xy_eval_upper = full(xy_eval_upper);
    xy_eval_lower = full(xy_eval_lower);
end
V_corridor = [xy_eval_upper,flip(xy_eval_lower,2),xy_eval_upper(:,1)]';

p1 = patch('Vertices',V_corridor,'Faces',1:length(V_corridor), ...
      'FaceColor','white','EdgeColor','red','FaceAlpha', 0.7);%,'Marker','.')%, 'LineStyle','--');

%% clean corridor
if 0
    inds_keep = closestPointsInCurve(xy_eval_upper);
    xy_eval_upper_clean = xy_eval_upper(:,inds_keep);
    inds_keep = closestPointsInCurve(xy_eval_lower);
    xy_eval_lower_clean = xy_eval_lower(:,inds_keep);
    
    V_corridor_clean = [xy_eval_upper_clean,flip(xy_eval_lower_clean,2),xy_eval_upper_clean(:,1)]';
    
    p2 = patch('Vertices',V_corridor_clean,'Faces',1:length(V_corridor_clean), ...
          'FaceColor','white','EdgeColor','black','FaceAlpha',0.7, ...
          'Linewidth',1,'LineStyle','--');%,'Marker','.','MarkerSize',0)%
end

end