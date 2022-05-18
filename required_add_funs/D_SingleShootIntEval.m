function WS_struct = D_SingleShootIntEval(OptsParams_WS,SetupI_WS,SetupII_WS)%(rhs,vars,integr_method,s,n,m,T,M)
    struct2CallerWS(OptsParams_WS)
    struct2CallerWS(SetupI_WS)
    struct2CallerWS(SetupII_WS)
    clearvars OptsParams_WS SetupI_WS SetupII_WS
    import casadi.*
    
    n_stages = s*n;

    % DT dynamics
    t    = SX.sym('t');
    t_f  = SX.sym('t_f');
    X_0  = SX.sym('X_0', n);
    U_0  = SX.sym('U_0', m*M);
    U_11 = SX.sym('U_11', m);
    w_L  = SX.sym('w_L', length(vars.C.w_L));
    Y    = [X_0; 0];
    h_vec_normed = SX.sym('h_vec_normed', M);

    Y_mat = [Y,repmat(Y,1,M)];
    Y_j   = Y;
    h_vec = h_vec_normed*t_f;
    U_0_mat = reshape(U_0,m,M);

    
    % preparation for different integration methods
    switch integr_method
        case {'IRK','SPRK'}
            STAGES = SX.sym('STAGES', n_stages*M);
            STAGES_mat = reshape(STAGES,n_stages,M);
        case 'IMP'
            X_interm = SX.sym('X_interm', n*M);
        case 'ERK'
            kappa     = SX.sym('kappa');
            rhs_Mayer = Function('rhs_Mayer', {vars.C.t,[vars.C.x;kappa],vars.C.u,vars.C.w_L}, {[rhs.Csym;L.Csym]});
    end
    
    % PATH FOLLOWING
    butcher_path = setButcherTable('ERK','4,1');
    Theta_0   = SX.sym('Theta_0', n_p);
    Theta_j   = Theta_0;
    Theta_mat = [Theta_0,repmat(Theta_0,1,M)];
    V_0       = SX.sym('V_0', m_p*M);
    V_0_mat   = reshape(V_0,m_p,M);


    %% integration on single shooting interval of length M
    for j=1:M
        X_0j     = Y_j(1:n);
        kappa_0j = Y_j(end);
        h = h_vec(j);

        % input
            U_0j = U_0_mat(:,j);
            if j<M
                U_0jp1 = U_0_mat(:,j+1);
            else
                U_0jp1 = U_11;
            end
            slope_U = (U_0jp1-U_0j)/h;
            if OC_U_ZOH==1
                slope_U = 0*slope_U;
            end
        
        switch integr_method
            case {'IRK','SPRK'}
                STAGES_j = STAGES_mat(:,j);%STAGES((j-1)*n_stages+1:j*n_stages);
                X_0jp1   = X_0j + h*weighted_stages_sum(STAGES_j);
                switch L_int_type
                    case 'Mayer'
                        kappa_0jp1 = kappa_0j + h*weighted_stages_kappa_sum(t,X_0j,U_0j,w_L,slope_U,h,STAGES_j);
                    case 'analytic'
                        kappa_0jp1 = kappa_0j + J_L_h_Casadi(X_0j,X_0jp1,U_0j,U_0jp1,h);
                    case 'GL'
                        STAGES_j_mat = reshape(STAGES_j,n,s);
%                       kappa_0jp1 = kappa_0j + exec_L_GL_int(n_GL,L_GL,h,X_0j,k_mat,butcher.b,U_0j,U_0jp1);
                        kappa_0jp1 = kappa_0j + exec_L_RK_time_int(L_GL,h,butcher.b,X_0j,STAGES_j_mat,U_0j,U_0jp1,w_L,n_GL);
                end
                Y_jp1 = [X_0jp1;kappa_0jp1];      
            case 'ERK'
                k_mat = compStagesRhsMat('explicit',[],rhs_Mayer,t,Y_j,U_0j,slope_U,h,butcher,{w_L});
                weighted_k_mat_sum = compWeightedStagesSum(k_mat,s,butcher.b);
                Y_jp1 = Y_j + h*weighted_k_mat_sum;

                switch L_int_type
                    case 'Mayer'
                        kappa_0jp1 = Y_jp1(end);
                    case 'analytic'
                        X_0jp1     = Y_jp1(1:end-1);
                        kappa_0jp1 = kappa_0j + J_L_h_Casadi(X_0j,X_0jp1,U_0j,U_0jp1,h);
                    case 'GL'
%                         kappa_0jp1 = kappa_0j + exec_L_GL_int(n_GL,L_GL,h,X_0j,k_mat(1:end-1,:),butcher.b,U_0j,U_0jp1);
                        kappa_0jp1 = kappa_0j + exec_L_RK_time_int(L_GL,h,butcher.b,X_0j,k_mat(1:end-1,:),U_0j,U_0jp1,w_L,n_GL);
                end
                Y_jp1(end) = kappa_0jp1;
            case 'IMP'
                X_0jp1     = X_interm((j-1)*n+1:j*n);                     
                k_x_1      = rhs.Cfun(t, 1/2*(X_0jp1+X_0j), U_0j+1/2*h*slope_U);
                kappa_0jp1 = kappa_0j + h*weighted_stages_kappa_sum(t,X_0j,U_0j,w_L,slope_U,h,k_x_1);

                Y_jp1 = [X_0j+h*k_x_1;
                         kappa_0jp1];
        end
        Y_mat(:,j+1) = Y_jp1;
        Y_j = Y_jp1;

        % PATH FOLLOWING (always use ERK)
        if PathFollowing==1
            k_mat_path = compStagesRhsMat('explicit',[],rhs_path.Cfun,t,Theta_j,V_0_mat(:,j),zeros(m_p,1),h,butcher_path,{});
            weighted_k_mat_sum_path = compWeightedStagesSum(k_mat_path,butcher_path.s,butcher_path.b);
            Theta_jp1 = Theta_j + h*weighted_k_mat_sum_path;
            Theta_mat(:,j+1) = Theta_jp1;
            Theta_j = Theta_jp1;
        end
    end


    X   = Y_jp1(1:end-1);
    J_L = Y_jp1(end);
    input_names   = {'t','X_0','U_0','U_11','h_vec_normed', 't_f','w_L'};
    input_values  = { t,  X_0,  U_0,  U_11,  h_vec_normed,   t_f,  w_L };
    output_names  = {'X_f','J_f'};
    output_values = { X,    J_L};
    switch integr_method
        case {'IRK','SPRK'}
            input_names   = horzcat(input_names(1:end-1), {'STAGES'},input_names(end));
            input_values  = horzcat(input_values(1:end-1),{ STAGES },input_values(end));
        case 'ERK'
            %do nothing
%             X_interm = Function('X_interm',{t,X_0,U_0,U_11,h_vec_normed,t_f},{X_interm});
        case 'IMP'
            input_names   = horzcat(input_names(1:end-1), {'X_interm'},input_names(end));
            input_values  = horzcat(input_values(1:end-1),{ X_interm },input_values(end));
    end

    F = Function('F',input_values,output_values, ...
                     input_names, output_names );
    X_mat = Function('X_mat',input_values(1:end-1),{Y_mat(1:n,:)});

    % PATH FOLLOWING
    if PathFollowing==1
        F_path    = Function('F_path',   {t,Theta_0,V_0,h_vec_normed,t_f},{Theta_jp1});
        Theta_mat = Function('Theta_mat',{t,Theta_0,V_0,h_vec_normed,t_f},{Theta_mat});
    end                           

    %%
    clearvars -except F t_f h N_sym X_interm tau_e h_vec X_mat ...
        F_path Theta_mat
    clearvars h
    strWS2WS_struct = who();
    strWS2WS_struct = strjoin(strWS2WS_struct,',');
    WS_struct = eval(sprintf('var2struct(%s);',strWS2WS_struct));
    try 
        WS_struct = rmfield(WS_struct,'ans');
    end
end