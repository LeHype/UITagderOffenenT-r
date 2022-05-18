function [xy_corridor,xy_corridor_plot,alpha_guess] = xyCorridorParam(xy_I1,t_grid_ed_eval,interp_method)

    syms theta [2 1] real
    theta_C = casadi.SX.sym('theta', 2);
    dxy_I1  = diff_pp(xy_I1);


    %% xy_corridor_plot
    xy_corridor_plot_smoothness = 'smooth';
%     xy_corridor_plot_smoothness = 'exact';

    xy_  = comp_pp_fun_sym( xy_I1,theta(1),xy_corridor_plot_smoothness);
    dxy_ = comp_pp_fun_sym(dxy_I1,theta(1),xy_corridor_plot_smoothness);
    
    alpha_ = angleBtwVec(dxy_,[1;0]);
    xy_corridor_ = xy_ + rot_z([0;theta(2)],-alpha_);
    xy_corridor_plot = matlabFunction(xy_corridor_,'Vars',{theta(1),theta(2)});

    alpha_guess_help = matlabFunction(alpha_,'Vars',{theta(1)});
    alpha_guess_fine = alpha_guess_help(t_grid_ed_eval);
    alpha_guess = griddedInterpolant(t_grid_ed_eval,alpha_guess_fine,interp_method);


    %% xy_corridor
    xy_  = comp_pp_fun_sym(xy_I1,theta(1),'smooth');
    dxy_ = comp_pp_fun_sym(dxy_I1,theta(1),'smooth');
    xy   = matlabFunction(xy_,'Vars',{theta(1)});
    dxy  = matlabFunction(dxy_,'Vars',{theta(1)});
    xy_C_   = xy(theta_C(1));
    dxy_C_  = dxy(theta_C(1));

    yo=1;
    switch yo
        case 1
            xy_pw  = pp2pwSymFun(xy_I1,theta(1));
            dxy_pw = pp2pwSymFun(diff_pp(xy_I1),theta(1));
        
            for ii=1:length(xy_pw.symfun)
                xy_ii_  =  xy_pw.symfun{ii};
                dxy_ii_ = dxy_pw.symfun{ii};
                xy_ii   = matlabFunction(xy_ii_,'Vars',{theta(1)});
                dxy_ii  = matlabFunction(dxy_ii_,'Vars',{theta(1)});
                xy_ii_C_  =  xy_ii(theta_C(1));
                dxy_ii_C_ = dxy_ii(theta_C(1));
                alpha_ii_C_ = angleBtwVec(dxy_ii_C_,[1;0],0);
                xy_corridor_ii_C_ = xy_ii_C_ + rot_z([0;theta_C(2)],-alpha_ii_C_);
                xy_corridor_pw_symfun_C(ii) = {xy_corridor_ii_C_};
            end
            xy_corridor_pw.symfun = xy_corridor_pw_symfun_C;
            xy_corridor_pw.intrvl = xy_pw.intrvl;
            xy_corridor_ = comp_p2_fun_sym(xy_corridor_pw,theta_C(1),'smooth');
        case 2
            alpha_C_ = angleBtwVec(dxy_C_,[1;0],0);
            xy_corridor_ = xy_C_ + rot_z([0;theta_C(2)],-alpha_C_);

%             alpha_ = angleBtwVec(dxy_,[1;0]);
%             xy_corridor_ = xy_ + rot_z([0;theta(2)],-alpha_);
    end

    xy_corridor = casadi.Function('xy_corridor',{theta_C},{xy_corridor_});
%     xy_corridor = matlabFunction(xy_corridor_,'Vars',{theta(1),theta(2)});



end