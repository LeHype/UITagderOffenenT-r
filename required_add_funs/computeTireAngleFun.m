function phi_fun = computeTireAngleFun(x_fun,tgrid)

TireRadius = 0.2;
x_traj = x_fun(tgrid);

interpolationMethod = 'spline';
% s_traj = [0, cumsum(abs(diff(vecnorm(x_traj(1:2,:)))))];
s_traj = -[0, cumsum(vecnorm( diff(x_traj(1:2,:),[],2) ))];
phi_traj = s_traj/TireRadius;
phi_fun_help = griddedInterpolant(tgrid,phi_traj',interpolationMethod,'previous');
phi_fun = @(t) phi_fun_help(t)';

end