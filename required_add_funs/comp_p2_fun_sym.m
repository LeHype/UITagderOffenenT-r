function pw_sym = comp_p2_fun_sym(pw,t,smoothness)
    
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
    
    
    N = length(pw.symfun);
    
    pw_sym = 0*t;
    for ii=1:N
        t_se   = pw.intrvl{ii};
        t_ii   = t_se(1);
        t_iip1 = t_se(2);
        pw_ii  = pw.symfun{ii};

        if ii==1
            pw_sym = pw_sym + pw_ii*(1 - HS(t-t_iip1));
        elseif ii>1 && ii<N
            pw_sym = pw_sym + pw_ii*SF(t_ii,t_iip1,t);
        elseif ii==N
            pw_sym = pw_sym + pw_ii*HS(t-t_ii);
        end
        
    end
end