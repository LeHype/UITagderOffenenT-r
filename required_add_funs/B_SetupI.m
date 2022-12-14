function WS_struct = B_SetupI(rhs,n,m,parametrize_OP,PathFollowing,QuantsOCP)
    arguments
        % mandatory
            rhs 
            n
            m
            parametrize_OP
            PathFollowing = 0
        % optional 
            % cost functional             
                QuantsOCP.L   = @(t,x,u,w_L)   w_L(1)*x'*eye(n)*x + w_L(2)*u'*eye(m)*u;
                QuantsOCP.Phi = @(t,x,u,w_Phi) 0;
                QuantsOCP.y   = @(t,x,u)       x;    
            % scale and weigh parts of cost functional   
                QuantsOCP.w_J           struct       = struct('L',ones(2,1),'Phi',0,'t_f',0)
                QuantsOCP.J_scale_coeff (1,3) double = [1,1,1]; %[c_J_t_f, c_J_L, c_J_Phi]   
            % t_f                                              
                QuantsOCP.lb_t_f    (1,1) double = 10^-3; % lower bound on final time t_f
                QuantsOCP.ub_t_f    (1,1) double = 100;   % upper bound on final time t_f
                QuantsOCP.t_f_guess (1,1) double = 10;  
            % choose T, h_des & M --> N                        
                QuantsOCP.T         (1,1) double = 10;
                QuantsOCP.h_des     (1,1) double = 10^-3; 
                QuantsOCP.M         {isinteger(QuantsOCP.M)} = 1;
                QuantsOCP.N         {isinteger(QuantsOCP.N)} = [];
            % x                                                
                % IC and general bounds 
                    QuantsOCP.x_0  (:,1) double =  ones(n,1);
                    QuantsOCP.lb_x (:,1) double = -Inf*ones(n,1);
                    QuantsOCP.ub_x (:,1) double =  Inf*ones(n,1);
                % constraint terminal state
                    QuantsOCP.constr_x_f        = 1; 
                    QuantsOCP.constr_states_x_f = 1:n;
                    QuantsOCP.box_states        = 1:n;
                    QuantsOCP.delta_box_x       = 0; %0 --> equality is forced
                    QuantsOCP.x_f               = 0*ones(n,1);
                % guess
                    QuantsOCP.x_guess_fun = [];
                    QuantsOCP.x_guess_fun_further_params = {};  
                % path constraints
                    QuantsOCP.x_path_constr = []; %Casadi Function of x_k
                    QuantsOCP.x_path_constr_lb = [];
                    QuantsOCP.x_path_constr_ub = [];
                % object avoidance
                    QuantsOCP.OA_constr    = []; %Casadi Function of x_k,OA_var 
                    QuantsOCP.OA_var_dim   = [];
                    QuantsOCP.OA_N_obj     = [];
                    % bound var
                        QuantsOCP.lb_OA_var    = [];
                        QuantsOCP.ub_OA_var    = [];
                    % bounds constraint
                        QuantsOCP.lb_OA_constr = []; %Casadi Function of x_k,OA_var 
                        QuantsOCP.ub_OA_constr = []; %Casadi Function of x_k,OA_var 
            % u                                                
                % general bounds
                    QuantsOCP.lb_u (:,1) double = -Inf*ones(m,1);
                    QuantsOCP.ub_u (:,1) double =  Inf*ones(m,1);
                % constraint initial and terminal input
                    QuantsOCP.constr_u_0  = 1;
                    QuantsOCP.constr_u_f  = 1;
                    QuantsOCP.delta_box_u = 0;
                % initial and terminal u
                    QuantsOCP.u_0 (:,1) double = 0*ones(m,1);
                    QuantsOCP.u_f (:,1) double = 0*ones(m,1);
                % guess
                    QuantsOCP.u_guess_fun = [];
                    QuantsOCP.u_guess_fun_further_params = {};        
            % du      
                % general bounds
                QuantsOCP.lb_du_s (:,1) double = -Inf*ones(m,1);
                QuantsOCP.ub_du_s (:,1) double =  Inf*ones(m,1);
                % guess
                QuantsOCP.du_guess_fun = [];
                QuantsOCP.du_guess_fun_further_params = {};            
            % y                                                
                QuantsOCP.lb_y = [];
                QuantsOCP.ub_y = [];                  
            % scaling of state and input                                        
                QuantsOCP.T_x_scale = [];
                QuantsOCP.T_u_scale = [];

            % PATH FOLLOWING
                % path-state dynamics
                    QuantsOCP.rhs_path = [];
                    QuantsOCP.n_p = [];
                    QuantsOCP.m_p = [];
                    QuantsOCP.xy_corridor = [];
                % theta 
                    %QuantsOCP.theta_0  = [xy_0 0 psi_0 0]';
                    QuantsOCP.lb_theta = [];
                    QuantsOCP.ub_theta = [];
                % v                                              
                    QuantsOCP.lb_v = [];
                    QuantsOCP.ub_v = [];
    end
    struct2CallerWS(QuantsOCP)
    import casadi.*

    n_w_L   = length(w_J.L);
    n_w_Phi = length(w_J.Phi);
%     t_f_guess = T;
    create_MatlabCasadi_sym('Matlab',{'t',sprintf('x %d 1',n),sprintf('u %d 1',m), ...
        sprintf('w_L %d 1',n_w_L),sprintf('w_Phi %d 1',n_w_Phi)},'caller')

    t_f_guess = T;
    
    %% dynamics and cost functionl                      
    rhs_ = rhs(t,x,u);
    L_   = L(t,x,u,w_L);
    Phi_ = Phi(t,x,u,w_Phi);
    y_   = y(t,x,u);
    
    %% x        
    x_f(setdiff(1:end,constr_states_x_f)) = x_0(setdiff(1:end,constr_states_x_f));
    x_0_num  = x_0;
    x_f_num  = x_f;
    if parametrize_OP==1
%         x_0_num  = x_0;
%         x_f_num  = x_f;
        x_0      = SX.sym('x_0',n);
        x_f      = SX.sym('x_f',n);
        x_0_orig = x_0;
        x_f_orig = x_f;
    end
    if isempty(x_guess_fun)
        x_guess_fun = @(t,x_0,x_f) x_0 + (x_f-x_0)/t_f_guess * t;
    end
    x_guess_norm_fun = @(tau,x_0,x_f) x_guess_fun(tau*t_f_guess,x_0,x_f);%,x_guess_fun_further_params{:}); % tau=t/t_f_guess <=> t=tau*t_f_guess, tau in [0;1]

    %% u   
    u_0_num = u_0;
    u_f_num = u_f;
    if parametrize_OP==1
%         u_0_num = u_0;
%         u_f_num = u_f;
        u_0 = SX.sym('u_0',m);
        u_f = SX.sym('u_f',m);
    end
    if isempty(u_guess_fun)
        u_guess_fun = @(t,u_0,u_f) u_0 + (u_f-u_0)/t_f_guess * t;
    end
    u_guess_norm_fun = @(tau,u_0,u_f) u_guess_fun(tau*t_f_guess,u_0,u_f);%,u_guess_fun_further_params{:}); % tau=t/t_f_guess <=> t=tau*t_f_guess, tau in [0;1]

    %% du
    if isempty(x_guess_fun)
        du_guess_fun = @(t,u_0,u_f) (u_f-u_0)/t_f_guess;
    end
    du_guess_norm_fun = @(tau,u_0,u_f) du_guess_fun(tau*t_f_guess,u_0,u_f); %,u_guess_fun_further_params{:}); % tau=t/t_f_guess <=> t=tau*t_f_guess, tau in [0;1]


    %% t_f
    t_f_in_CostFun = 1*(w_J.t_f>0);
    if parametrize_OP==1
        lb_t_f_num = lb_t_f;
        ub_t_f_num = ub_t_f;
        lb_t_f     = SX.sym('lb_t_f',1);
        ub_t_f     = SX.sym('ub_t_f',1);
        t_f_guess_num = t_f_guess;
        t_f_guess     = SX.sym('t_f_guess',1);
    end

    %% T, h_des & M --> N                        
    if parametrize_OP==1
        import casadi.*
    end  
    if isempty(N)
        N = ceil(1/M * T/(h_des)); %number of control intervals --> determine N s.t. N*M*h = T <=> h = T/(N*M) (h<=h_des is as close as possible to h_des)
    else
        fprintf(2,[repmat('_',1,60),'\n'])
        warning('N was chosen. This means that T, h_des & M are not involved in determining h. Hence the stepsize h depends on your choice of N solely.')
        fprintf(2,[repmat('¯',1,60),'\n'])
    end
    
    %% determine w_J_infos
    w_J_fns = fieldnames(w_J);
    for ii=1:length(w_J_fns)
        w_J_dims(ii) = length(w_J.(w_J_fns{ii}));
    end
    w_J_inds_help = cumsum(w_J_dims);
    for ii=1:length(w_J_fns)%-1
        if ii==1
            w_J_inds(ii) = {1:w_J_inds_help(ii)};
        else
            w_J_inds(ii) = {(w_J_inds_help(ii-1)+1) : w_J_inds_help(ii)};
        end
    end
    w_J_infos = [w_J_fns';w_J_inds;struct2cell(w_J)'];
    clearvars w_J_fns w_J_dims W_J_inds_help w_J_inds

    %% PATH FOLLOWING
    if PathFollowing==0
        n_p = 1;
        m_p = 1;
        rhs_path = @(t,theta,v) 0*theta;
    end
    create_MatlabCasadi_sym('Matlab',{sprintf('theta %d 1',n_p),sprintf('v %d 1',m_p)},'caller')
    rhs_path_ = rhs_path(t,theta,v);

    %% output WS_struct
    clearvars rhs parametrize_OP QuantsOCP t x u L Phi y ans 
    clearvars rhs_path theta v
    strWS2WS_struct = who();
    strWS2WS_struct = strjoin(strWS2WS_struct,',');
    WS_struct = eval(sprintf('var2struct(%s);',strWS2WS_struct));
    try 
        WS_struct = rmfield(WS_struct,'ans');
    end

end