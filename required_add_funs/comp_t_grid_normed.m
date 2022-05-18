function t_grid_normed = comp_t_grid_normed(t_grid_normed_opt,N,M)

    switch t_grid_normed_opt
        case 1
            t_grid_normed = linspace(0,1,N*M+1);
        case 2
            t_grid_normed = linspace(0,10,N*M+1).^2 / 10^2;
        case 3
            use_h_min = 0;
            h_min = 2*1/N*M;
            h_max = 2.5*1/N*M;
            if h_min<(1/N*M) || h_max<(1/N*M)
                error('Choice of h_min cannot be fulfilled')
            elseif h_max==1
                error('Choice of h_max cannot be fulfilled')
            elseif h_max>=0.3 
                warning('Choice of h_max is quite large')
            end
            d = @(alpha) (1-exp(-alpha))/(N*M);
            if use_h_min==1
                f_alpha = @(alpha,h_min) -1/alpha*log(1-d(alpha)) - 0 - h_min;
                b_h = h_min;
            else
                f_alpha = @(alpha,h_max) 1 - (-1/alpha*log(exp(-alpha)+d(alpha))) - h_max;
                b_h = h_max;
            end
            alpha_num = fzero(@(alpha) f_alpha(alpha,b_h),1);
            y_grid = linspace(exp(-alpha_num),1,N*M+1);
            fun = @(y) -1/alpha_num*log(y);  
            t_grid_normed = flip(fun(y_grid));
    end

end