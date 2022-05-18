function WS_struct = G_PostProcessOCP2(sol,X_mat,BuildOCP_WS,OptsParams_WS,n,m,M,N)%,xy_corridor)
    arguments
    % mandatory
        sol
        X_mat
        BuildOCP_WS struct
        OptsParams_WS struct
        n (1,1) double
        m (1,1) double
        M (1,1) double
        N (1,1) double
%         xy_corridor
    % optional
        % nothing
    end
    struct2CallerWS(OptsParams_WS)
    struct2CallerWS(BuildOCP_WS)
    import casadi.*

    % extract data
        w_opt     = full(sol.x);
        lam_w_opt = full(sol.lam_x);
        lam_g_opt = full(sol.lam_g);
        J_opt     = full(sol.f);


    t_f_opt    = w_opt(end);
    w_parts    = get_w_parts(w_opt,n,m,M,~OC_U_ZOH);
    w_mid      = w_parts.middle;
    w_mid_rs   = reshape(w_mid,[],N);
    U_traj_OCP = w_mid_rs(1:M*m,:);
    u_traj_OCP = reshape(U_traj_OCP,m,[]);
    x_traj_OCP = w_mid_rs(M*m+1:M*m+n,:);

    Theta_traj = w_mid_rs(end-3:end-2,:);
    V_traj     = w_mid_rs(end-1:end,:);

    if strcmp(integr_method,'ERK')==0
        R_traj_OCP = w_mid_rs(M*m+n+1:end,:);
        r_traj_OCP = reshape(R_traj_OCP,[],M);
    end

    
%     if strcmp(integr_method,'ERK_adaptive')==1
%         H_traj = w_mid_rs(m+n+1:end,:);
%         h_traj = reshape(H_traj,1,[]);
%         tgrid_u = cumsum(h_traj);
%         tgrid_x = tgrid_u; %doesnt work for M>1
%     else
        tgrid_u = t_grid_normed*t_f_opt;
        tgrid_x = tgrid_u;
        help = reshape(tgrid_u(1:end-1)',M,N);
        tgrid_x_OCP = help(1,:); 
%     end
    
    if OC_U_ZOH==1
        u_traj = [u_traj_OCP, u_traj_OCP(:,end)];
        plot_u = @stairs;
        u_fun = griddedInterpolant(tgrid_u,u_traj','previous');
    else
        u_traj = [u_traj_OCP, w_parts.extra];
        plot_u = @plot;
        u_fun = griddedInterpolant(tgrid_u,u_traj','linear');
    end
    
    h_mat_normed = reshape(diff(t_grid_normed),M,N);

    if M>1
        x_traj = [];
        for ii=1:size(x_traj_OCP,2)
            if ii<size(x_traj_OCP,2)
                U_11 = U_traj_OCP(1:m,ii+1);
                X_mat_cols = 2:M;
            else
                U_11 = u_traj(:,end);
                X_mat_cols = 2:M+1;
            end
            if strcmp(integr_method,'ERK')==1
                X_mat_ii = X_mat(tgrid_x_OCP(ii),x_traj_OCP(:,ii),U_traj_OCP(:,ii),U_11,h_mat_normed(:,ii),t_f_opt);
            else
                X_mat_ii = X_mat(tgrid_x_OCP(ii),x_traj_OCP(:,ii),U_traj_OCP(:,ii),U_11,h_mat_normed(:,ii),t_f_opt,R_traj_OCP(:,ii));
            end
            
            x_traj = [x_traj, x_traj_OCP(:,ii), X_mat_ii(:,X_mat_cols)];
        end
    else
        U_11 = u_traj(:,end);
        ii = size(x_traj_OCP,2);
        if strcmp(integr_method,'ERK')==1
            X_mat_ii = X_mat(tgrid_x_OCP(ii),x_traj_OCP(:,ii),U_traj_OCP(:,ii),U_11,h_mat_normed(:,ii),t_f_opt);
        else
            X_mat_ii = X_mat(tgrid_x_OCP(ii),x_traj_OCP(:,ii),U_traj_OCP(:,ii),U_11,h_mat_normed(:,ii),t_f_opt,R_traj_OCP(:,ii));
        end
        x_traj = [x_traj_OCP,X_mat_ii(:,end)];
    end

    x_traj = full(x_traj);
    x_fun = griddedInterpolant(tgrid_x,x_traj','linear');


%     MakeDefaultFig(100,100,'Screen',1)
%     subplot(1,2,1)
%         plot(tgrid_x(1:end-1),Theta_traj,'Marker','.')
%         legend('$\theta_1$','$\theta_2$')
%     subplot(1,2,2)
%         plot(tgrid_x(1:end-1),V_traj)
%         legend('$v_1$','$v_2$')

    %%
%     clearvars  
    clearvars(strjoin(fieldnames(OptsParams_WS),' ')) 
    clearvars(strjoin(fieldnames(BuildOCP_WS),' ')) 
    clearvars OptsParams_WS BuildOCP_WS

    strWS2WS_struct = who();
    strWS2WS_struct = strjoin(strWS2WS_struct,',');
    WS_struct = eval(sprintf('var2struct(%s);',strWS2WS_struct));
    try 
        WS_struct = rmfield(WS_struct,'ans');
    end

end