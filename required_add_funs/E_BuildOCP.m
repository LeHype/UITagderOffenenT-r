function WS_struct = E_BuildOCP(SingleShootIntEval_WS,OptsParams_WS,SetupI_WS,SetupII_WS)
    struct2CallerWS(SingleShootIntEval_WS)
    struct2CallerWS(OptsParams_WS)
    struct2CallerWS(SetupI_WS)
    struct2CallerWS(SetupII_WS)
    clearvars SingleShootIntEval_WS OptsParams_WS SetupI_WS SetupII_WS
    import casadi.*

%     constr_X_interm = 1
    
    %% init
    % Start with an empty NLP
        J_L = 0;
        structStructureFieldnames = {'U','X','X_interm','X_path','dU','STAGES','h_vec','Theta','V','OA'};
        var_name_cell = {'w','w_0','lb_w','ub_w','g','lb_g','ub_g'};        
        % create struct "structStructure" with fields "structStructureFieldnames" that are empty
            help = repmat({[]},size(structStructureFieldnames));
            help = reshape([structStructureFieldnames;help],1,[]);
            structStructure = struct(help{:});
        % create struct "var_name_struct" with fields "var_name_cell" that are initialize with structStructure
            help = repmat({structStructure},size(var_name_cell));
            help = reshape([var_name_cell;help],1,[]);
            var_name_struct = struct(help{:});
        % unpack struct "var_name_struct"
            struct2CallerWS(var_name_struct)
            % --> w,w_0,lb_w,ub_w,g,lb_g,ub_g with fields U,X,dU,STAGES (that are empty)     
    % initial conditions
        X_0 = SX.sym('X_0', n);
        U_0 = SX.sym('U_0', m*M);
    % PATH FOLLOWING
        if PathFollowing==1
            Theta_0 = SX.sym('Theta_0', n_p);
            V_0     = SX.sym('V_0', m_p*M);
        end
    % obstacle avoidance
        if ~isempty(OA_constr) && ~isempty(OA_var_dim)
            OA_var_0 = SX.sym('OA_var_0', OA_N_obj*OA_var_dim);
        end
    % parameter vector
        % w_J
            N_p_w_J = w_J_infos{2,end};
            p_w_J = SX.sym('p_w_J', N_p_w_J);
            for ii=1:size(w_J_infos,2)
                p_w_J_parts.(['w_',w_J_infos{1,ii}]) = p_w_J(w_J_infos{2,ii}); 
            end
            struct2CallerWS(p_w_J_parts);
            p_w_J_num = cell2mat(w_J_infos(3,:)');
        % h_vec
            p_h_mat  = SX.sym('p_h_mat',M,N);
            p_h_vecs = reshape(p_h_mat,[M*N,1]);
%             t_grid_normed_opt = 1
            t_grid_normed = comp_t_grid_normed(t_grid_normed_opt,N,M);
            p_h_vecs_num = diff(t_grid_normed)';
        % stack parts of p
            p     = [p_w_J; p_h_vecs];
            p_num = [p_w_J_num; p_h_vecs_num];


    %% initial guess
    x_guess_norm_fun_t  = @(tau) x_guess_norm_fun(tau,x_0_num,x_f_num);
    u_guess_norm_fun_t  = @(tau) u_guess_norm_fun(tau,u_0_num,u_f_num);
    du_guess_norm_fun_t = @(tau) u_guess_norm_fun(tau,u_0_num,u_f_num);

    %% U,X: initial constraint
    x_eq = zeros(n,1);

    X_start = SX.sym('X_start', n);
    U_start = SX.sym('U_start', M*m);
    w_first    = [U_start; X_start];
    lb_w_first = [[zeros((M-1)*m,1);lb_u_0]; lb_x_0];
    ub_w_first = [[zeros((M-1)*m,1);ub_u_0]; ub_x_0];
    w_0_first  = [[zeros((M-1)*m,1);u_0]; x_0];

    X_first_left  = X_start;
    U_first_left  = U_start((M-1)*m+1:M*m);
    U_first_right = U_0(1:m);
    X_first_right = X_0;
    g_first    = [U_first_left-U_first_right; X_first_left-X_first_right];
    lb_g_first = [zeros(m,1); x_eq];
    ub_g_first = [zeros(m,1); x_eq];
    

    w.X    = [w.X,    X_0]; %free initial u(0) for reform_tilde_x=1 && constr_u_0==1
    lb_w.X = [lb_w.X, lb_x];
    ub_w.X = [ub_w.X, ub_x];
    w_0.X  = [w_0.X,  x_0];

    % PATH FOLLOWING
        if PathFollowing==1
            w.Theta    = [w.Theta,    Theta_0];
            lb_w.Theta = [lb_w.Theta, lb_theta];
            ub_w.Theta = [ub_w.Theta, ub_theta];
%             lb_w.Theta = [lb_w.Theta, [-5*10^-2;lb_theta(2)]]; % due to arctan Approximation
%             ub_w.Theta = [ub_w.Theta, [ 5*10^-2;ub_theta(2)]];
            w_0.Theta  = [w_0.Theta,  zeros(n_p,1)];

%             g.Theta    = [g.Theta,    xy_corridor(Theta_0)-X_0([1,2])];
%             lb_g.Theta = [lb_g.Theta, zeros(n_p,1)];
%             ub_g.Theta = [ub_g.Theta, zeros(n_p,1)];

            J_path = norm(xy_corridor(Theta_0)-X_0([1,2]));
        end
    % obstacle avoidance 
        if ~isempty(OA_constr) && ~isempty(OA_var_dim)
            w.OA    = [w.OA,    OA_var_0];
            lb_w.OA = [lb_w.OA, repmat(lb_OA_var,OA_N_obj,1)];
            ub_w.OA = [ub_w.OA, repmat(ub_OA_var,OA_N_obj,1)];
            w_0.OA  = [w_0.OA,  repmat(lb_OA_var + (ub_OA_var-ub_OA_var)/2,OA_N_obj,1)];

            g.OA    = [g.OA,    OA_constr(X_0,OA_var_0)];
            lb_g.OA = [lb_g.OA, lb_OA_constr];
            ub_g.OA = [ub_g.OA, ub_OA_constr];
        end

   

    %% +++ for loop: U, X (and STAGES) +++
    warning_X_k = 0;
    for k=0:N-1
        %% k+1 -> k
        if k>0
            X_k   = X_kp1_right;
            U_km1 = U_k;
            U_k   = U_kp1;
            if PathFollowing==1
                Theta_k = Theta_kp1_right;
                V_k     = V_kp1;
            end
            if ~isempty(OA_constr) && ~isempty(OA_var_dim)
                OA_var_k = OA_var_kp1;
            end
        else
            X_k   = X_0;
            U_km1 = U_first_left; % --> U_first_left-U_first_right = U_first_left-U_{0,1} != 0 
                                  % --> U_km1 = U_{0,1}
                                  % --> diff_U = [U_{0,1}-U_{0,1}; U_{0,1}-U_{0,2}; ...] 
                                  %            = [0;               U_{0,1}-U_{0,2}; ...]
                                  % --> first element always satisfies b_du_s-constraint
            U_k = U_0;
            if PathFollowing==1
                Theta_k = Theta_0;
                V_k     = V_0;
            end
            if ~isempty(OA_constr) && ~isempty(OA_var_dim)
                OA_var_k = OA_var_0;
            end
        end

        %% create variables
        % U: next control bundle
            if k~=N-1
                U_kp1 = SX.sym(['U_' num2str(k+1)], m*M);
            else
                U_kp1 = SX.sym('U_last_left', m);
            end
            if OC_U_ZOH==1
                U_kp11 = zeros(m,1); 
            else
                U_kp11 = U_kp1(1:m); 
            end
        % X: next state at end of shooting interval
            X_kp1_right = SX.sym(['X_' num2str(k+1)], n);                
        % STAGES/X_interm/h_vec
            h_vec_normed_k = p_h_mat(:,k+1);
            h_vec_k = h_vec_normed_k*t_f;
            switch integr_method
                case {'IRK','SPRK'}
                    STAGES_k = SX.sym(['STAGES_' num2str(k)], n_stages*M);     
                case 'ERK'
                    % do nothing
                case 'IMP'
                    if M>1
                        X_interm_k_help = SX.sym(['X_interm_',num2str(k)],n*(M-1));
                        X_interm_k = [X_interm_k_help; X_kp1_right];    
                    else
                        X_interm_k = X_kp1_right;
                    end
            end
        % PATH FOLLOWING
            if PathFollowing==1
                Theta_kp1_right = SX.sym(['Theta_' num2str(k+1)], n_p);
                V_kp1 = SX.sym(['V_' num2str(k+1)], m_p*M);
            end
        % obstacle avoidance 
            if ~isempty(OA_constr) && ~isempty(OA_var_dim)
                OA_var_kp1 = SX.sym(['OA_var_' num2str(k+1)], OA_N_obj*OA_var_dim);
            end


        %% Integrate till the end of the interval
        if k>0
            t_now = sum(sum(p_h_mat(:,1:k)));
        else
            t_now = 0;
        end
        NVPs = {'t',t_now,  'X_0',X_k, 'U_0',U_k, 'U_11',U_kp11, 'w_L',w_L, 'h_vec_normed',h_vec_normed_k, 't_f',t_f};
        switch integr_method
            case {'IRK','SPRK'}
                NVPs = horzcat(NVPs, {'STAGES',STAGES_k}); 
            case 'ERK'
                % do nothing                
            case 'IMP'
                NVPs = horzcat(NVPs, {'X_interm',X_interm_k});  
        end
        F_k = F(NVPs{:});
        J_k_end    = F_k.J_f;
        X_kp1_left = F_k.X_f;
        J_L        = J_L + J_k_end;  


        if PathFollowing==1
            Theta_kp1_left = F_path(t_now,Theta_k,V_k,h_vec_normed_k,t_f);
            if M>1
                Theta_k_mat = Theta_mat(t_now,Theta_k,V_k,h_vec_normed_k,t_f);
                Theta_k_interm_mat = Theta_k_mat(:,2:end-1);
            end
        end



        %% impose constraints on variables
        if 1
            %% X
            if k<N-1
                w.X    = [w.X,    X_kp1_right];
                lb_w.X = [lb_w.X, lb_x];
                ub_w.X = [ub_w.X, ub_x];
    
                % Add equality constraint: dynamics
                g.X    = [g.X,    X_kp1_left-X_kp1_right]; % g(x)==0 <=>
                lb_g.X = [lb_g.X, x_eq];                   % g(x)>=0 AND 
                ub_g.X = [ub_g.X, x_eq];                   % g(x)<=0

                if ~isempty(x_path_constr)
                    g.X_path    = [g.X_path,    x_path_constr(X_kp1_right)];
                    lb_g.X_path = [lb_g.X_path, x_path_constr_lb];
                    ub_g.X_path = [ub_g.X_path, x_path_constr_ub];
                end
            end

            %% U
            w.U  = [w.U, U_k];
            if k==N-1 && OC_U_ZOH==1
                lb_w.U = [lb_w.U, [repmat(lb_u,M-1,1); lb_u_f]];
                ub_w.U = [ub_w.U, [repmat(ub_u,M-1,1); ub_u_f]];
            else
                lb_w.U = [lb_w.U, repmat(lb_u,M,1)];
                ub_w.U = [ub_w.U, repmat(ub_u,M,1)];
            end    

            
            %% delta_U
            if reform_dyn_xu==0 && any(abs([lb_du_s,ub_du_s])~=Inf,'all')
                diff_U = diff(reshape([U_km1(end-m+1:end);U_k],m,M+1),1,2); %=diff(U_{k-1,M},U_{k,1},...,U_{k,M})
                %      lb_du_s*h <= diff_U   <= ub_du_s*h
                % <=>  lb_du_s   <= diff_U/h <= ub_du_s
                help = reshape(diff_U'./h_vec_k,m*M,1);
                g.dU    = [g.dU, help];
                lb_g.dU = [lb_g.dU, repmat(lb_du_s,M,1)];
                ub_g.dU = [ub_g.dU, repmat(ub_du_s,M,1)];
            end       

            %% PATH FOLLOWING: Theta and V
            if PathFollowing==1 %&& k<N-1
                if k<N-1
                    w.Theta    = [w.Theta,     Theta_kp1_right];
                    lb_w.Theta = [lb_w.Theta,  lb_theta];
                    ub_w.Theta = [ub_w.Theta,  ub_theta];

                    % Add equality constraint: dynamics
                    g.Theta    = [g.Theta,    Theta_kp1_left-Theta_kp1_right]; % g(x)==0 <=>
                    lb_g.Theta = [lb_g.Theta, zeros(n_p,1)];                   % g(x)>=0 AND 
                    ub_g.Theta = [ub_g.Theta, zeros(n_p,1)];                   % g(x)<=0

%                     J_path = J_path + norm(xy_corridor(Theta_kp1_right)-X_kp1_right([1,2]));
                end
                if k<N-1 % && k>0 
                    % Add 'quasi'-equality constraint: stay inside corrido
                    g.Theta    = [g.Theta,    xy_corridor(Theta_kp1_right)-X_kp1_right([1,2])];
                    lb_g.Theta = [lb_g.Theta, zeros(n_p,1)];
                    ub_g.Theta = [ub_g.Theta, zeros(n_p,1)];   
%                     lb_g.Theta = [lb_g.Theta, -5*10^-3*ones(n_p,1)];
%                     ub_g.Theta = [ub_g.Theta,  5*10^-3*ones(n_p,1)]; 
                end


                w.V    = [w.V, V_k];
                lb_w.V = [lb_w.V, repmat(lb_v,M,1)];
                ub_w.V = [ub_w.V, repmat(ub_v,M,1)];
            end

            %% obstacle avoidance
            if k>0
                if ~isempty(OA_constr) && ~isempty(OA_var_dim)
                    w.OA    = [w.OA,    OA_var_kp1];
                    lb_w.OA = [lb_w.OA, repmat(lb_OA_var,OA_N_obj,1)];
                    ub_w.OA = [ub_w.OA, repmat(ub_OA_var,OA_N_obj,1)];
                    w_0.OA  = [w_0.OA,  repmat(lb_OA_var + (ub_OA_var-ub_OA_var)/2,OA_N_obj,1)];
    
                    g.OA    = [g.OA,    OA_constr(X_k,OA_var_kp1)];
                    lb_g.OA = [lb_g.OA, lb_OA_constr];
                    ub_g.OA = [ub_g.OA, ub_OA_constr];
                end
            end


            %% STAGES/X_interm/h_vec
            switch integr_method
                case {'IRK','SPRK'}
                    w.STAGES    = [w.STAGES,     STAGES_k];
                    lb_w.STAGES = [lb_w.STAGES, -Inf*ones(n_stages*M,1)];
                    ub_w.STAGES = [ub_w.STAGES,  Inf*ones(n_stages*M,1)];

                    X_kj = X_k;
                    g_help    = [];
                    lb_g_help = [];
                    ub_g_help = [];
                    for j=1:M
                        h = h_vec_k(j);
                        [U_kj, slope_U] = determine__U_kj__slope_U(U_k,U_kp1,j,k,m,M,N,h,OC_U_ZOH);                        
                        STAGES_kj = STAGES_k((j-1)*n_stages+1:j*n_stages); 

                        % state fixpoint eq. 
                            stages_FP_kj = stages_FP(STAGES_kj,0,X_kj,U_kj,slope_U,h);                            
                            g_help    = [g_help;    stages_FP_kj];
                            lb_g_help = [lb_g_help; stages_eq];
                            ub_g_help = [ub_g_help; stages_eq];

                        X_kj = X_kj + h*weighted_stages_sum(STAGES_kj);
                    end
                    g.STAGES    = [g.STAGES,    g_help];
                    lb_g.STAGES = [lb_g.STAGES, lb_g_help];
                    ub_g.STAGES = [ub_g.STAGES, ub_g_help];

                    if warning_X_k~=1
                        warning('Did not constraint intermediate X_k yet')
                    end
                    warning_X_k = 1;
                case 'IMP'
                    if M>1
                        w.STAGES    = [w.STAGES,    X_interm_k_help];
                        lb_w.STAGES = [lb_w.STAGES, repmat(lb_x,M-1,1)];
                        ub_w.STAGES = [ub_w.STAGES, repmat(ub_x,M-1,1)];
                        
                        X_kj = X_k;
                        g_help    = [];
                        lb_g_help = [];
                        ub_g_help = [];
                        for j=1:M-1
                            [U_kj, slope_U] = determine__U_kj__slope_U(U_k,U_kp1,j,k,m,M,N,h,OC_U_ZOH);                        
                            X_kjp1 = X_interm_k_help((j-1)*n+1:j*n); 
                            k_x_1  = rhs.Cfun(t_now,1/2*(X_kjp1+X_kj),U_kj + 1/2*h*slope_U);
                            X_kjp1_rhs = X_kj + h*k_x_1;

                            % state fixpoint eq. 
                                g_help    = [g_help;    X_kjp1_rhs-X_kjp1];
                                lb_g_help = [lb_g_help; x_eq];
                                ub_g_help = [ub_g_help; x_eq];

                            % j+1 -> j
                            X_kj = X_kjp1_rhs;
                        end
                        g.STAGES    = [g.STAGES,    g_help];
                        lb_g.STAGES = [lb_g.STAGES, lb_g_help];
                        ub_g.STAGES = [ub_g.STAGES, ub_g_help];                         
                    end                       
                case 'ERK'
                    if constr_X_interm==1 && M>1
                        X_mat_k = X_mat(t_now,X_k,U_k,U_kp11,h_vec_normed_k,t_f);
                        g.X_interm    = [g.X_interm,    reshape(X_mat_k(:,2:end-1),numel(X_mat_k(:,2:end-1)),1)]; %X_interm(0,X_k,U_k,U_kp11,h_vec_normed,t_f)];
                        lb_g.X_interm = [lb_g.X_interm, repmat(lb_x,M-1,1)];
                        ub_g.X_interm = [ub_g.X_interm, repmat(ub_x,M-1,1)];   

                        if ~isempty(x_path_constr)
                            for ii=2:M-1
                                g.X_path    = [g.X_path,    x_path_constr(X_mat_k(:,ii))];
                                lb_g.X_path = [lb_g.X_path, x_path_constr_lb];
                                ub_g.X_path = [ub_g.X_path, x_path_constr_ub];
                            end
                        end

                    end
            end
            

            %% 
        end

        %% initial guess of variables
        if 1
            t_guess_1 = ((k+1)*M)/(N*M);     %normalized time                 at k+1
            t_guess_2 = (k*M+(0:M-1))/(N*M); %normalized time vector starting at k

            %% X
            if k<N-1
                w_0.X  = [w_0.X,  x_guess_norm_fun_t(t_guess_1)];
            end
            %% U 
            if reform_dyn_xu==0
                w_0.U  = [w_0.U, reshape(u_guess_norm_fun_t(t_guess_2),[],1)];  
            else
                w_0.U  = [w_0.U, reshape(du_guess_norm_fun_t(t_guess_2),[],1)];
            end

            %% PATH FOLLOWING
            if PathFollowing==1 
                if k<N-1
                    w_0.Theta = [w_0.Theta, k/N*[1;0]];
                end
                w_0.V = [w_0.V, zeros(m_p*M,1)];
            end

                      
            %% STAGES/X_interm/h_vec
            switch integr_method
                case {'IRK','SPRK'}
                    w_0.STAGES  = [w_0.STAGES, zeros(n_stages*M,1)];
                case 'IMP'
                    if M>1
                        w_0.STAGES  = [w_0.STAGES,  reshape(x_guess_norm_fun_t(t_guess_2(2:end)),[],1)];                      
                    end                       
                case 'ERK'
                    % do nothing
            end
        end
    end
    X_last_left = X_kp1_left;


    %% U: U_k at t_k=t_f for OC_U_ZOH==0
    if OC_U_ZOH==0
        %U
            fprintf('Add extra U_k at t=t_f\n')
            U_last_left = U_kp1(1:m);
            w_extra    = U_last_left;
            lb_w_extra = lb_u;
            ub_w_extra = ub_u;
            if reform_dyn_xu==0
                w_0_extra = reshape(u_guess_norm_fun_t(1),[],1);  
            else
                w_0_extra = reshape(du_guess_norm_fun_t(1),[],1);
            end
        %dU
            if reform_dyn_xu==0 && any(abs([lb_du_s,ub_du_s])~=Inf,'all')
                diff_U  = diff(reshape([U_k(end-m+1:end);U_last_left],m,2),1,2); %=diff(U_{k,M},U_{k+1,1},...,U_{k+1,M})
                h = h_vec_k(end);
                help    = reshape(diff_U/h,m,1);
                g_extra    = help;
                lb_g_extra = lb_du_s;
                ub_g_extra = ub_du_s; 
            else
                g_extra    = [];
                lb_g_extra = [];
                ub_g_extra = [];
            end
    else
        U_last_left = U_k(end-m+1:end);
        w_extra    = [];
        lb_w_extra = [];
        ub_w_extra = [];
        w_0_extra  = [];  
        g_extra    = [];
        lb_g_extra = [];
        ub_g_extra = [];
    end

    %% U,X: terminal constraint
    X_end = SX.sym('X_end', n);
    U_end = SX.sym('U_end', m);
    w_last    = [U_end;  X_end];
    lb_w_last = [lb_u_f; lb_x_f];
    ub_w_last = [ub_u_f; ub_x_f];
    w_0_last  = [u_f;    x_f];
    
    X_last_right = X_end;
    U_last_right = U_end;
%     g.X    = [g.X,    X_last_left-X_last_right]; % here "k==N-1", i.e. last iter. of FL
%     lb_g.X = [lb_g.X, x_eq];                    
%     ub_g.X = [ub_g.X, x_eq];   

    g_last    = [U_last_left-U_last_right; X_last_left-X_last_right];
    lb_g_last = [zeros(m,1); x_eq];
    ub_g_last = [zeros(m,1); x_eq];


    %% concatenate contents cells
    n_var_U      = numel(w.U);
    n_var_X      = numel(w.X);
    n_var_STAGES = numel(w.STAGES);
    n_con_X      = numel(g.X);
    n_con_dU     = numel(g.dU);
    n_con_STAGES = numel(g.STAGES);
    for ii=1:length(var_name_cell)
%         var_name_cell{ii}
        evalc(['cur_var = ',var_name_cell{ii},';']);
        if vars_as_block == 1 %--> w = [U,X,STAGES]'
            for mm=1:length(structStructureFieldnames)
                cur_field = structStructureFieldnames{mm};
                if ~isempty(cur_var.(cur_field))
                    cur_var.(cur_field) = reshape(cur_var.(cur_field),numel(cur_var.(cur_field)),1);
                else
                    cur_var = rmfield(cur_var,cur_field);
                end
            end
            cur_var = struct2cell(cur_var);
            cur_var = vertcat(cur_var{:});
        else
            cur_var_help = [];
            max_dim2 = 0;
            for rr=1:length(structStructureFieldnames)
                max_dim2 = max([max_dim2, size(cur_var.(structStructureFieldnames{rr}),2)]);
            end
            for ll=1:max_dim2
                for mm=1:length(structStructureFieldnames)
                    cur_field = structStructureFieldnames{mm};
                    try
                        cur_var_help = [cur_var_help; cur_var.(cur_field)(:,ll)];
                    catch
                        if strcmp(cur_field,'OA')==1
                            hi=1;
                        end
                    end
                end
            end
            cur_var = cur_var_help;
        end
        evalc([var_name_cell{ii}, ' = cur_var;']);            
    end

    %% add final time t_f to variables to be optimized
    if t_f_in_CostFun==1
        lb_t_f_if  = lb_t_f;
        ub_t_f_if  = ub_t_f;
        w_0_t_f_if = t_f_guess;
    else
        lb_t_f_if  = T;
        ub_t_f_if  = T;
        w_0_t_f_if = T;
    end
    w    = [w_first;    w;    w_extra;    w_last;    t_f       ];
    lb_w = [lb_w_first; lb_w; lb_w_extra; lb_w_last; lb_t_f_if ];
    ub_w = [ub_w_first; ub_w; ub_w_extra; ub_w_last; ub_t_f_if ];
    w_0  = [w_0_first;  w_0;  w_0_extra;  w_0_last;  w_0_t_f_if];
    g    = [g_first;    g;    g_extra;    g_last   ];      
    lb_g = [lb_g_first; lb_g; lb_g_extra; lb_g_last]; 
    ub_g = [ub_g_first; ub_g; ub_g_extra; ub_g_last];


    %% cost functional
    J_L_fun = Function('J_L_fun',{w,w_L},{J_L});
    J_t     = w_t_f*t_f;
%     J_Phi   = Phi.Cfun(N_sym*h,X_last_right,U_last_right,w_Phi);
        J_Phi   = 0*t_f;%Phi.Cfun(t_f,X_last_right,U_last_right,w_Phi);
    if t_f_in_CostFun==1
        J_int      =   J_scale_coeff(1) * J_L;
        J_terminal =   J_scale_coeff(2) * J_Phi ...
                     + J_scale_coeff(3) * J_t;
    else
        J_int      = J_scale_coeff(1) * J_L;
        J_terminal = J_scale_coeff(2) * J_Phi;
    end
    J = J_int + J_terminal;
    if PathFollowing==1 
        J = J + 10^-6*J_path;
        fprintf('Added path costs\n\n');
    end

%     t_prob_form_f = datevec(datetime('now'));
%     prob_form_time = etime(t_prob_form_f,t_prob_form_0);
%     fprintf('Required time to formulate problem: %0.2f sec.\n',prob_form_time )

    fprintf('Number of optim. variables: %d   ---   Number of constraints: %d\n',length(w),length(g))


    %% parametrized OP
    if parametrize_OP==1
        pOP_text_cell = {'w_0','lb_w','ub_w','lb_g','ub_g'};
        for ii=1:length(pOP_text_cell)
            cur_var = pOP_text_cell{ii};
            sprintf_help = repmat({cur_var},1,3);
            sym_vars_help = '{x_0_orig,x_f_orig,u_0,u_f,t_f_guess,lb_t_f,ub_t_f}';
            num_vars_help = '(x_0_num,x_f_num,u_0_num,u_f_num,t_f_guess_num,lb_t_f_num,ub_t_f_num)';
            evalc(sprintf( ['%s_fun = Function(''%s_fun'', ',sym_vars_help,',{%s});'], ...
                           sprintf_help{:}) );
            evalc(sprintf( ['%s     = %s_fun',num_vars_help,';'], ...
                           sprintf_help{1:2}) );
        end
    end


    %%
    clearvars -except J J_int J_terminal w_0 w lb_w ub_w g lb_g ub_g p ...
        n_var_U n_var_X n_var_STAGES n_con_X n_con_dU n_con_STAGES ...
        N M diff_tau_e_all_fun tau_e_all_fun ...
        t_grid_normed_opt t_grid_normed p p_w_J p_h_vecs p_num p_w_J_num p_h_vecs_num
    strWS2WS_struct = who();
    strWS2WS_struct = strjoin(strWS2WS_struct,',');
    WS_struct = eval(sprintf('var2struct(%s);',strWS2WS_struct));
    try 
        WS_struct = rmfield(WS_struct,'ans');
    end
end