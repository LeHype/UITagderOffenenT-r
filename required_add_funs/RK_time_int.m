function F = RK_time_int(f,param,f_0,t_0,t_f,h,n_GL)%,b,c)

%%
if 0
    N = round((t_f-t_0)/h);
    s = length(b);
    F = f_0;
    if ~isempty(param)
        k_mat = repmat(0*f(0,param{:}),1,s);
    else
        k_mat = repmat(0*f(0),1,s);
    end
    
    for k=0:N-1
        t_k = t_0 + k*h;
        for ii=1:s
            if ~isempty(param)
                k_mat(:,ii) = f(t_k+c(ii)*h,param{:});
            else
                k_mat(:,ii) = f(t_k+c(ii)*h);
            end
        end
        weighted_k_mat     = (diag(b)*k_mat')';
        weighted_k_mat_sum = sum(weighted_k_mat,2);
        F = F + h*weighted_k_mat_sum;
    end
end

%%
[b_GL,c_GL] = get_GL_Int_coeffs(n_GL); %coeffecients for time integral from -1 to 1

N = round((t_f-t_0)/h);
s = length(b_GL);
F = f_0;
if ~isempty(param)
    k_mat = repmat(0*f(0,param{:}),1,s);
else
    k_mat = repmat(0*f(0),1,s);
end

t = linspace(t_0,t_f,N+1);

for k=1:N
    a = t(k);
    b = t(k+1);    
    t_b_GL = (b-a)/2*b_GL;
    for ii=1:s
        t_eval_i = (b-a)/2*c_GL(ii) + (a+b)/2; %integral transformation s.t. it goes from a to b
        if ~isempty(param)
            k_mat(:,ii) = f(t_eval_i,param{:});
        else
            k_mat(:,ii) = f(t_eval_i);
        end
    end
    weighted_k_mat     = (diag(t_b_GL)*k_mat')';
    weighted_k_mat_sum = sum(weighted_k_mat,2);
    F = F + weighted_k_mat_sum;
end


end