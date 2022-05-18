function pw = pp2pwSymFun(pp,t)  

    N = pp.pieces;
    dim = pp.dim;
    for ii=1:N
        t_ii   = pp.breaks(ii);
        t_iip1 = pp.breaks(ii+1);
        coefs_ii = pp.coefs((ii-1)*dim+1:ii*dim,:);
        
        for jj=1:dim
            pp_ii_help(jj,1) = poly2sym(coefs_ii(jj,:),t);
        end
        pp_ii_ = subs(pp_ii_help, t, t-t_ii);
        pw_symfun(ii) = {pp_ii_};
        pw_intrvl(ii) = {[t_ii,t_iip1]};
    end
    pw.symfun = pw_symfun;
    pw.intrvl = pw_intrvl;
end