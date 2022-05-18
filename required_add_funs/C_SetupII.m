function WS_struct = C_SetupII(OptsParams_WS,SetupI_WS)
    struct2CallerWS(OptsParams_WS)
    struct2CallerWS(SetupI_WS)
    clearvars OptsParams_WS SetupI_WS
    keep_vars = '';
    
    %%
    [butcher,s] = setButcherTable(integr_method,des_int_order,integr_method_params);

    if deltaEMP_Mayer_on==1
        [butcher_kappa,s] = setButcherTable('ERK','2,3',integr_method_params);
        warning('Uses DeltaEMP method')
    end
    if ~exist('butcher_kappa')
        butcher_kappa = butcher;
    end

%     fprintf('\n')
    keep_vars = [keep_vars,var2varnameString(butcher,s)];

    %% COMPUTE rhs, L, Phi, y and vars
    OCP_rhs_L_data = set_OCP_sysdims_vars_scaling(n,m,n_w_L,n_w_Phi,scale_x,scale_u,T_x_scale,T_u_scale,PathFollowing,n_p,m_p);
    OCP_rhs_L_data.rhs_ = rhs_;
    OCP_rhs_L_data.L_   = L_;
    OCP_rhs_L_data.L_   = L_;
    OCP_rhs_L_data.Phi_ = Phi_;
    OCP_rhs_L_data.y_   = y_;
    OCP_rhs_L_data.rhs_path_ = rhs_path_;
    [rhs,L,Phi,y,vars,n,m,rhs_path] = compute_rhs_L(OCP_rhs_L_data,reform_dyn_xu,PathFollowing); % n=n+m for reform_dyn_xu==1
    % rhs and L have inputs t,x,u
    keep_vars = [keep_vars,var2varnameString(rhs,L,Phi,y,vars,n,m,rhs_path,n_p,m_p)];
    
    %% x
    lb_x_0 = x_0;
    ub_x_0 = x_0;
%     x_f_constr = x_f;%(constr_states_x_f);
    if constr_x_f==0
        lb_x_f = -Inf*ones(n,1);%x_f_constr;
        ub_x_f =  Inf*ones(n,1);%x_f_constr;
    else
        lb_x_f = x_f;
        ub_x_f = x_f;

        box_x_f_help1 = eye(length(x_f(box_states)));
        box_x_f_help2 = diag(sign(x_f(box_states)))*delta_box_x;
        lb_x_f(box_states) = (box_x_f_help1-box_x_f_help2) * x_f(box_states);
        ub_x_f(box_states) = (box_x_f_help1+box_x_f_help2) * x_f(box_states);


        lb_x_f(box_states) = x_f(box_states) - delta_box_x;
        ub_x_f(box_states) = x_f(box_states) + delta_box_x;
    end
    unconstr_states = setdiff(1:n,constr_states_x_f);
    if ~isempty(unconstr_states)
        lb_x_f(unconstr_states) = -Inf*ones(length(unconstr_states),1);
        ub_x_f(unconstr_states) =  Inf*ones(length(unconstr_states),1);
    end
    x_eq = zeros(n,1);
    keep_vars = [keep_vars,var2varnameString(lb_x_0,ub_x_0,lb_x_f,ub_x_f,x_eq)];


    %% u
    if constr_u_0==1
        box_u_0_help = diag(sign(u_0))*delta_box_u;
        lb_u_0 = (eye(m)-box_u_0_help) * u_0;
        ub_u_0 = (eye(m)-box_u_0_help) * u_0;
    else
        lb_u_0 = lb_u;
        ub_u_0 = ub_u;   
    end
    if constr_u_f==1
        box_u_f_help = diag(sign(u_f))*delta_box_u;
        lb_u_f = (eye(m)-box_u_f_help) * u_f;
        ub_u_f = (eye(m)-box_u_f_help) * u_f; 
    else
        lb_u_f = lb_u;
        ub_u_f = ub_u;   
    end
    keep_vars = [keep_vars,var2varnameString(lb_u_0,ub_u_0,lb_x_f,ub_x_f)];

    %% stages
    n_stages  = s*n;
    stages_eq = zeros(n_stages,1);
    keep_vars = [keep_vars,var2varnameString(n_stages,stages_eq)];


    %% reform_dyn_xu
    if reform_dyn_xu==1
        x_0    = [x_0;u_0];
        lb_x   = [lb_x;lb_u];
        ub_x   = [ub_x;ub_u];
        lb_x_0 = [lb_x_0;lb_u_0];
        ub_x_0 = [ub_x_0;ub_u_0];
        lb_x_f = [lb_x_f;lb_u_f];
        ub_x_f = [ub_x_f;ub_u_f];
        x_guess_norm_fun = @(tau,x_0,x_f) [x_guess_norm_fun(tau,x_0(1:n-m),x_f(1:n-m)); u_guess_norm_fun(tau,x_0(n-m+1:n),x_f(n-m+1:n))];
        constr_states_x_f = [constr_states_x_f,n];

        lb_u   = lb_du_s;
        ub_u   = ub_du_s;
        lb_u_0 = lb_du_s;
        ub_u_0 = ub_du_s;
        lb_u_f = lb_du_s;
        ub_u_f = ub_du_s;
        u_guess_norm_fun = du_guess_norm_fun;

        lb_du_s = -Inf; % Delta(Delta U) is not constraint
        ub_du_s =  Inf; % Delta(Delta U) is not constraint
        
        % update system order and equality vectors
        n_stages  = s*n;
        x_eq      = zeros(n,1);
        stages_eq = zeros(n_stages,1);
    end
    keep_vars = [keep_vars,var2varnameString(x_0,lb_x,ub_x,lb_x_0,ub_x_0,lb_x_f,ub_x_f,x_guess_norm_fun,constr_states_x_f, ...
                    lb_u,ub_u,lb_u_0,ub_u_0,lb_u_f,ub_u_f,u_guess_norm_fun, ...
                    lb_du_s,ub_du_s,n_stages,x_eq,stages_eq)];


    %% compute_stages_FP
    switch integr_method
%         case 'SPRK'%,'IMP'}
%             [stages_FP,weighted_stages_sum,stages,weighted_stages_kappa_sum] = ...
%                 compute_stages_FP_SPRK_Casadi(butcher,rhs_SPRK_,t,x_SPRK,u,L.Csym); %[rhs_SPRK_;L_energy_]
%             % --> stages_FP(stages,t,x,u,h)
        case 'ERK'
            % do nothing
        case {'IRK','IMP'}
            warning('Why is here butcher_kappa?')
            [stages_FP,weighted_stages_sum,stages,weighted_stages_kappa_sum] = ...
                compute_stages_FP_Casadi(butcher_kappa,rhs.Csym,vars.C,L.Csym);
%             [stages_FP,weighted_stages_sum,stages,weighted_stages_kappa_sum] = ...
%                 compute_stages_FP_Casadi(butcher,rhs.Csym,vars.C,L.Csym);
    end
    if ismember(integr_method,{'SPRK','IRK','IMP'})
        keep_vars = [keep_vars,var2varnameString(stages_FP,weighted_stages_sum,stages,weighted_stages_kappa_sum)];
    end

    %% try to compute analytic solution to cost integral (x and u are linearly interpolated)
    switch L_int_type
        case 'Mayer'
            % do nothing
        case 'analytic'
            J_L_h_Casadi = compute_analyt_in_cost(L.Mfun,n,m);
%             % test
%                 kappa_add = J_L_h_Casadi('x_0',0.5*ones(n,1),'x_1',0.8*ones(n,1),'u_0',1000,'u_1',1500,'h',0.1);
%                 disp(full(kappa_add.res))
%                 if isnan(full(kappa_add.res))
%                     warning('NaN in cost function! analyt_int_cost is set to zero!')
%                     L_int_type = 'Mayer';
%                 end
        case 'GL'
            L_GL = create_L_GL(L.Cfun,n,m,n_w_L);
    end
    keep_vars = [keep_vars,'J_L_h_Casadi L_GL '];

    %%
%     tic
    evalc(['clearvars -except ',keep_vars]);
    strWS2WS_struct = who();
    strWS2WS_struct = strjoin(strWS2WS_struct,',');
    WS_struct = eval(sprintf('var2struct(%s);',strWS2WS_struct));
    try 
        WS_struct = rmfield(WS_struct,'ans');
    end
%     toc

end