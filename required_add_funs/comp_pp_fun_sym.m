function pp_fun_sym = comp_pp_fun_sym(pp_fun,t,smoothness)
    
    if nargin<3
        smoothness = 'exact';
    end
    
    switch smoothness
        case 'exact'
            HS = @(t)     heaviside(t);
            SF = @(t,a,b) stepfunction(t,a,b);
        case 'smooth'
            HS = @(t)     smoothHeaviside(t);
            SF = @(t,a,b) smoothStepFunction(t,a,b);
    end
    
    
    N = pp_fun.pieces;
    dim = pp_fun.dim;
    
    pp_fun_sym = 0*t;
    for ii=1:N
        t_ii   = pp_fun.breaks(ii);
        t_iip1 = pp_fun.breaks(ii+1);
        
        coefs_ii = pp_fun.coefs((ii-1)*dim+1:ii*dim,:);
        
        for jj=1:dim
            poly_ii_help(jj,1) = poly2sym(coefs_ii(jj,:),t);
        end
        poly_ii = subs(poly_ii_help, t, t-t_ii);
%         if ii>1
            pp_fun_sym = pp_fun_sym + poly_ii*SF(t_ii,t_iip1,t);
%         else
%             pp_fun_sym = pp_fun_sym + poly_ii*(1 - HS(t-t_iip1));
%         end
    end
    pp_fun_sym = pp_fun_sym + ppval(pp_fun,t_iip1)*HS(t-t_iip1); %pp_fun_sym(t)=pp_fun(t_f) for all t>=t_f
end