function J_L_h_Casadi = compute_analyt_in_cost(L,n,m)
    n_str = int2str(n);
    m_str = int2str(m);
    create_MatlabCasadi_sym('Matlab',{'t','h',['x_0 ',n_str,' 1'],['x_1 ',n_str,' 1'],['u_0 ',m_str,' 1'],['u_1 ',m_str,' 1']},'caller')

    x_slope = (x_1-x_0)/h;
    u_slope = (u_1-u_0)/h;

    x_traj = x_0 + x_slope*t;
    u_traj = u_0 + u_slope*t;

    if 1
        L_FOH_       = L(t,x_traj,u_traj);
        J_L_h_indef_ = expand(int(L_FOH_,t,'IgnoreAnalyticConstraints',true));
        J_L_h_indef  = matlabFunction(J_L_h_indef_,'Vars',{t,x_0,x_1,u_0,u_1,h});
        J_L_h_       = J_L_h_indef(h,x_0,x_1,u_0,u_1,h); %energy cost at t=0 is zero %- J_energy_h_indef(0,x_0,x_1,u_0,u_1,h)       
%         J_L_h_       = J_L_h_indef(h,x_0,x_1,u_0,u_1,h) - J_L_h_indef(0,x_0,x_1,u_0,u_1,h)
    else
        J_L_h_ = int(L(t,x_traj,u_traj),t,0,h,'IgnoreAnalyticConstraints',true)
    end
    % done because int 1/x does not give log(|x|) but log(x) and x can be negative
    % by (Matlabs) def. log(x) = log(|x|) + e^(i*x) and therefore real(log(x)) = log(|x|) 
        J_L_h_ = real(J_L_h_);
%         z = log(abs(x_0));
% %         z = abs(x_0);
%         J_L_h_ = transpose(z)*z
    J_L_h  = matlabFunction(J_L_h_,'Vars',{x_0,x_1,u_0,u_1,h});

    create_MatlabCasadi_sym('Casadi',{'t','h',['x_0 ',n_str,' 1'],['x_1 ',n_str,' 1'],['u_0 ',m_str,' 1'],['u_1 ',m_str,' 1']},'caller')
    J_L_h_Casadi_ = J_L_h(x_0,x_1,u_0,u_1,h);
    import casadi.*
%     J_L_h_Casadi  = Function('J_L_h',{x_0,x_1,u_0,u_1,h},{J_L_h_Casadi_},{'x_0','x_1','u_0','u_1','h'},{'res'});
    J_L_h_Casadi  = Function('J_L_h',{x_0,x_1,u_0,u_1,h},{J_L_h_Casadi_});


end