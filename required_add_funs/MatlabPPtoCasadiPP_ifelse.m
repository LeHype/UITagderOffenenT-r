function [pp_C,pp_C_,theta] = MatlabPPtoCasadiPP_ifelse(pp_M)
% pp_M is object obtained via interp1(t,x,interp_method,'pp')
% https://github.com/casadi/casadi/issues/2489

% !!! derivative of MatlabPPtoCasadiPP (different approach than if_else) evaluated at breaks yields nan !!! :(


import casadi.*

dim            = pp_M.dim;
order          = pp_M.order;
breaks_vec     = pp_M.breaks;
coeffs_mat_all = pp_M.coefs;

theta = SX.sym('theta',1,1);

theta_shift_help = theta-breaks_vec(1:end-1);
theta_exps       = flip(0:order-1);

pp_C_ = SX.zeros(dim,1);
% % keep value at start/end constant
% for ii=1:dim
%     pp_C_(ii,1) = if_else(theta<breaks_vec(1), coeffs_mat_all(ii,:)*(0.^theta_exps)', 0);
%     pp_C_(ii,1) = pp_C_(ii,1) + ...
%                   if_else(theta>=breaks_vec(end), coeffs_mat_all(end-dim+ii,:)*(diff(breaks_vec(end-1:end)).^theta_exps)', 0);
% end
for ii=1:dim
    coeffs_mat = coeffs_mat_all(ii:dim:size(coeffs_mat_all,1),:);
    for jj=1:length(breaks_vec)-1
        pp_jj = coeffs_mat(jj,:)*(theta_shift_help(jj).^theta_exps)';
        if jj==1
            %extrapoltaion to -Inf with pp at the start
            cond = theta<breaks_vec(jj+1);
        elseif jj==length(breaks_vec)-1 
            %extrapoltaion to +Inf with pp at the end
            cond = theta>=breaks_vec(jj);
        else
            cond = theta>=breaks_vec(jj) & theta<breaks_vec(jj+1);
        end
        pp_C_(ii,1) = pp_C_(ii,1) + if_else(cond,pp_jj,0);
    end
end

pp_C = Function('pp_C',{theta},{pp_C_});

end