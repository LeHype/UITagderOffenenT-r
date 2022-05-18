function [stages_FP,weighted_stages_sum,stages,weighted_stages_kappa_sum] = ...
            compute_stages_FP_Casadi(butcher,rhs_,vars,L_)

import casadi.*

%% Butcher
b = butcher.b;
c = butcher.c;
s = length(c);


%% prepare sym. variables for stages
t   = vars.t;
x   = vars.x;
u   = vars.u;
w_L = vars.w_L;
n = length(rhs_);
m = length(u);

rhs = Function('rhs',{t,x,u},    {rhs_});
L   = Function('L',  {t,x,u,w_L},{L_});

stages  = SX.sym('stages',n*s,1);
u_slope = SX.sym('u_slope',m,1);
h       = SX.sym('h');

stages_mat = reshape(stages,n,s);


%% stages_FP
stages_rhs_mat = compStagesRhsMat('implicit',stages_mat,rhs,t,x,u,u_slope,h,butcher,{});
stages_rhs = reshape(stages_rhs_mat,numel(stages_rhs_mat),1);
stages_FP_ = stages_rhs - stages; % [k_1_rhs;...,k_s_rhs] - [k_1;...;k_s] == 0
stages_FP  = Function('stages_FP',{stages,t,x,u,u_slope,h},{stages_FP_});


%% weighted_stages_sum & weighted_stages_kappa_sum
weighted_stages_sum_ = compWeightedStagesSum(stages_mat,s,b);
weighted_stages_sum  = Function('weighted_stages_sum',{stages},{weighted_stages_sum_});

stages_kappa_mat           = compStagesRhsMat('implicit',stages_mat,L,t,x,u,u_slope,h,butcher,{w_L});
weighted_stages_kappa_sum_ = compWeightedStagesSum(stages_kappa_mat,s,b);
weighted_stages_kappa_sum  = Function('weighted_stages_kappa_sum',{t,x,u,w_L,u_slope,h,stages}, {weighted_stages_kappa_sum_});

end



