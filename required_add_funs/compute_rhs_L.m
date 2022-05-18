function [rhs_OP,L_OP,Phi_OP,y_OP,vars,n,m,rhs_path_OP] = compute_rhs_L(OCP_rhs_L_data,reform_dyn_xu,PathFollowing)

%% allocate variables
struct2CallerWS(OCP_rhs_L_data)

if isa(L_,'sym')==0
    L_ = sym(0);
end
if isa(Phi_,'sym')==0
    Phi_ = sym(0);
end

%% scaling
if 1
    %% rhs
    % unscaled
        rhs  = matlabFunction(rhs_,'Vars',{t,x,u});
    % scaled
        % d/dt x = f(t,x,u)  
        %   <==> T_x * d/dt x = T_x * f(t,T_x^-1*T_x*x,T_u^-1*T_u*u)   z:=T_x*x, v:=T_u*u
        %   <==> d/dt z = T_x*f(t,T_x^-1*z,T_u^-1*v)
        rhs_ = T_x_scale * rhs(t,inv(T_x_scale)*x,inv(T_u_scale)*u);
        rhs  = matlabFunction(rhs_,'Vars',{t,x,u});

    %% running cost L
    % unscaled
        L  = matlabFunction(L_,'Vars',{t,x,u,w_L});
    % scaled
        L_ = L(t,inv(T_x_scale)*x,inv(T_u_scale)*u,w_L);
        if isa(L_,'sym')==0
            L_ = sym(0);
        end
        L  = matlabFunction(L_,'Vars',{t,x,u,w_L});

    %% terminal cost Phi
    % unscaled
        Phi  = matlabFunction(Phi_,'Vars',{t,x,u,w_Phi});
    % scaled
        Phi_ = Phi(t,inv(T_x_scale)*x,u,w_Phi);
        if isa(Phi_,'sym')==0
            Phi_ = sym(0);
        end
        Phi  = matlabFunction(Phi_,'Vars',{t,x,u,w_Phi});

    %% output y
    % unscaled
        y  = matlabFunction(y_,'Vars',{t,x,u});
    % scaled
        y_ = y(t,inv(T_x_scale)*x,inv(T_u_scale)*u);% + 0*t; %0*t to keep it symbolic if
        y  = matlabFunction(y_,'Vars',{t,x,u});

    %% PATH FOLLOWING
        rhs_path = matlabFunction(rhs_path_,'Vars',{t,theta,v});
        
end

%% reform_dyn_xu
rhs_Mfun_orig = rhs;
L_Mfun_orig   = L;
Phi_Mfun_orig = Phi;
y_Mfun_orig   = y;
if reform_dyn_xu==1   
    create_MatlabCasadi_sym('Matlab',{sprintf('x %d 1',n+m),sprintf('u %d 1',m)},'caller')
    % x_new = [x_old;u_old], u_new = delta_u_old
    
    rhs_ = [rhs(t,x(1:n),x(n+1:n+m)); u];
    rhs  = matlabFunction(rhs_,'Vars',{t,x,u});

    L_ = L(t,x(1:n),x(n+1:n+m),w_L);
    if isa(L_,'sym')==0
        L_ = sym(0);
    end
    L  = matlabFunction(L_,'Vars',{t,x,u,w_L}); % new u=du does not appear but is kept due to consistency
    
    Phi_ = Phi(t,x(1:n),x(n+1:n+m),w_Phi);
    if isa(Phi_,'sym')==0
        Phi_ = sym(0);
    end
    Phi  = matlabFunction(Phi_,'Vars',{t,x,u,w_Phi}); % new u=du does not appear but is kept due to consistency
    
    y_ = y(t,x(1:n),x(n+1:n+m));
    y  = matlabFunction(y_,'Vars',{t,x,u}); % new u=du does not appear but is kept due to consistency
    
    %overwrite n
    n = n+m;
end

%% variables
% Matlab
    vars.M.t     = t;
    vars.M.x     = x;
    vars.M.u     = u;
    vars.M.w_L   = w_L;
    vars.M.w_Phi = w_Phi;
    vars.M.y     = y_;
    vars.M.theta = theta;
    vars.M.v     = v;
% Casadi
    create_MatlabCasadi_sym('Casadi',{'t',sprintf('x %d 1',n),sprintf('u %d 1',m), ...
        sprintf('w_L %d 1',length(w_L)),sprintf('w_Phi %d 1',length(w_Phi)), ...
        sprintf('theta %d 1',n_p),sprintf('v %d 1',m_p)},'caller')
    y_ = y(t,x,u);
    vars.C.t     = t;
    vars.C.x     = x;
    vars.C.u     = u;
    vars.C.w_L   = w_L;
    vars.C.w_Phi = w_Phi;
    vars.C.y     = y_;
    vars.C.theta = theta;
    vars.C.v     = v;
    
%% Casadi formatting
import casadi.*
% rhs
    rhs_Casadi_ = rhs(t,x,u);
    rhs_Casadi  = Function('rhs_Casadi',{t,x,u},{rhs_Casadi_});
% running cost L
    L_Casadi_ = L(t,x,u,w_L);
    L_Casadi  = Function('L_Casadi',{t,x,u,w_L},{L_Casadi_});
% terminal cost Phi
    Phi_Casadi_ = Phi(t,x,u,w_Phi);
    Phi_Casadi  = Function('Phi_Casadi',{t,x,u,w_Phi},{Phi_Casadi_});
% output y
    y_Casadi_ = y(t,x,u);
    y_Casadi  = Function('y_Casadi',{t,x,u},{y_Casadi_});
% PATH FOLLOWING
    if PathFollowing==1
        rhs_path_Casadi_ = rhs_path(t,theta,v);    
        rhs_path_Casadi  = Function('rhs_path_Casadi',{t,theta,v},{rhs_path_Casadi_});
    else
        rhs_path_Casadi_ = [];
        rhs_path_Casadi  = [];
    end

    
%% pack results in struct
% rhs
    rhs_OP.Mfun = rhs;
    rhs_OP.Msym = rhs_;
    rhs_OP.Cfun = rhs_Casadi;
    rhs_OP.Csym = rhs_Casadi_;
% running cost L
    L_OP.Mfun = L;
    L_OP.Msym = L_;
    L_OP.Cfun = L_Casadi;
    L_OP.Csym = L_Casadi_;
% terminal cost Phi
    Phi_OP.Mfun = Phi;
    Phi_OP.Msym = Phi_;
    Phi_OP.Cfun = Phi_Casadi;
    Phi_OP.Csym = Phi_Casadi_;
% terminal cost Phi
    y_OP.Mfun = y;
    y_OP.Msym = y_;
    y_OP.Cfun = y_Casadi;
    y_OP.Csym = y_Casadi_;

% if reform_dyn_xu==1
    rhs_OP.Mfun_orig = rhs_Mfun_orig;
    L_OP.Mfun_orig   = L_Mfun_orig;
    Phi_OP.Mfun_orig = Phi_Mfun_orig;
    y_OP.Mfun_orig   = y_Mfun_orig;

% PATH FOLLOWING
    rhs_path_OP.Mfun = rhs_path;
    rhs_path_OP.Msym = rhs_path_;
    rhs_path_OP.Cfun = rhs_path_Casadi;
    rhs_path_OP.Csym = rhs_path_Casadi_;

end

