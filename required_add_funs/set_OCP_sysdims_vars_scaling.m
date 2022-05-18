function OCP_rhs_L_data = set_OCP_sysdims_vars_scaling(n,m,n_w_L,n_w_Phi,scale_x,scale_u,T_x_scale,T_u_scale, ...
                                                            PathFollowing,n_p,m_p)
    create_MatlabCasadi_sym('Matlab',{'t',sprintf('x %d 1',n),sprintf('u %d 1',m), ...
        sprintf('w_L %d 1',n_w_L),sprintf('w_Phi %d 1',n_w_Phi)},'caller')
    
    % scaling matrices
    if isempty(T_x_scale)==1 || scale_x==0
        T_x_scale = eye(n);
    end
    if isempty(T_u_scale)==1 || scale_u==0
        T_u_scale = eye(m);
    end

%     if PathFollowing==0
%         OCP_rhs_L_data = var2struct(n,m,t,x,u,w_L,w_Phi,T_x_scale,T_u_scale);
%     else
        create_MatlabCasadi_sym('Matlab',{sprintf('theta %d 1',n_p),sprintf('v %d 1',m_p)},'caller')
        OCP_rhs_L_data = var2struct(n,m,t,x,u,w_L,w_Phi,T_x_scale,T_u_scale,n_p,m_p,theta,v);
%     end

end
