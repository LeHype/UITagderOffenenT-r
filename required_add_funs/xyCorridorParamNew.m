function [xy_corridor,xy_corridor_plot,alpha_guess] = xyCorridorParamNew(xy_I1,t_grid_ed_eval,interp_method)
    import casadi.*

    [pp_C,pp_C_,v] = MatlabPPtoCasadiPP_ifelse(xy_I1);
    dpp_C_ = jacobian(pp_C_,v);
%     dpp_C = Function('dpp_C',{v},{dpp_C_});

    
    theta = casadi.SX.sym('theta',1,1);
    
    alpha_ = angleBtwVec(dpp_C_,[1;0],0);
    pp_C_corridor_ = pp_C_ + rot_z([0;theta],-alpha_);
    xy_corridor  = Function('xy_corridor',{[v;theta]},{pp_C_corridor_});

    xy_corridor_plot = xy_corridor;

    alpha_guess_help = Function('alpha_guess',{v},{alpha_});
    alpha_guess_fine = full(alpha_guess_help(t_grid_ed_eval));
    alpha_guess = griddedInterpolant(t_grid_ed_eval,alpha_guess_fine,interp_method);

end