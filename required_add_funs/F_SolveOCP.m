function WS_struct = F_SolveOCP(BuildOCP_WS,OptsParams_WS,n,m,options)
    arguments
    % mandatory
        BuildOCP_WS   struct
        OptsParams_WS struct
        n             (1,1) double
        m             (1,1) double
    % optional
%         options.casadi_options struct = [];
        options.AdjointInitialGuess {mustBeMember(options.AdjointInitialGuess,[0,1])} = 0;
        options.lam_w_0 = [];
        options.lam_g_0 = [];
    end
    struct2CallerWS(OptsParams_WS)
    struct2CallerWS(BuildOCP_WS)
    struct2CallerWS(options)
    import casadi.*


    %% solve the OCP
    % Create an NLP solver %TODO QP, QCQP, ...
        t_prep_0 = datevec(datetime('now'));
        prob = struct('f', J, 'x', w, 'g', g, 'p', p);
        if ~isempty(casadi_options)
            solver = nlpsol('solver', solver_choice, prob, casadi_options);
        else
            solver = nlpsol('solver', solver_choice, prob);
        end
        t_prep_f = datevec(datetime('now'));
        prep_time = etime(t_prep_f,t_prep_0);
    % set numerical values (bounds and initial guess
        solver_fun_NVPs = {'x0',  w_0, ...
                           'p',   p_num, ...
                           'lbx', lb_w, 'ubx', ub_w, ...
                           'lbg', lb_g, 'ubg', ub_g};
        if AdjointInitialGuess==1 && (~isempty(lam_g_0) && ~isempty(lam_w_0))
            solver_fun_NVPs = [solver_fun_NVPs, ...
                               'lam_g0', lam_g_0, ...
                               'lam_x0', lam_w_0];
        elseif AdjointInitialGuess==1 && (isempty(lam_g_0) || isempty(lam_w_0))
            error('Supply both lam_g_0 and lam_w_0.')
        end
    % solve the OCP
        t_comp_0 = datevec(datetime('now'));
        sol = solver(solver_fun_NVPs{:});
        t_comp_f = datevec(datetime('now'));
        sol_time = etime(t_comp_f,t_comp_0);
    % extract data
        w_opt     = full(sol.x);
        lam_w_opt = full(sol.lam_x);
        lam_g_opt = full(sol.lam_g);
        J_opt     = full(sol.f);


    %%
%     evalc(['clearvars -except ',keep_vars]);
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