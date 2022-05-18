function J_L_all = exec_L_RK_time_int(L,h,b,x_k,k_mat,u_k,u_kp1,w_L,n_GL)
    s = length(b);
    x_kip1_mat = repmat(0*x_k,1,s);
    u_kip1_mat = repmat(0*u_k,1,s);
    du_k       = (u_kp1-u_k)/h;
    
    x_ki = x_k;
    u_ki = u_k;
    for ii=1:s-1
        %compute
            x_kip1 = x_ki + h*b(ii)*k_mat(:,ii);
            u_kip1 = u_ki + h*b(ii)*du_k;
        %save
            x_kip1_mat(:,ii) = x_kip1;
            u_kip1_mat(:,ii) = u_kip1;
        %update
            x_ki = x_kip1;
            u_ki = u_kip1;
    end
    dx_ki_mat = k_mat;
    
    x_ki_mat = [x_k, x_kip1_mat];
    u_ki_mat = [u_k, u_kip1_mat]; 
    
    import casadi.*
    create_MatlabCasadi_sym('Casadi',{'t'},'caller','SX')

    J_L_all = 0;
    for ii=1:s
        param = {h,x_ki_mat(:,ii),dx_ki_mat(:,ii),u_ki_mat(:,ii),du_k,w_L};        
        J_L_all = J_L_all + RK_time_int(L,param,0,  0,  b(ii),b(ii),n_GL);%,b_L,c_L);
                           %RK_time_int(f,param,f_0,t_0,t_f,  h,         b,  c  );
    end
    J_L_all = h*J_L_all;    
end

