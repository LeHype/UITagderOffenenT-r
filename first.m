function [] = first()

warning('off','all');
% commandwindow
%             import java.awt.*;
%             import java.awt.event.*;
%             %Create a Robot-object to do the key-pressing
%             rob=Robot;
%              rob.keyPress(KeyEvent.VK_CONTROL);
%               rob.keyPress(KeyEvent.VK_SHIFT);
%                rob.keyPress(KeyEvent.VK_G)
%              rob.keyRelease(KeyEvent.VK_CONTROL)
%            rob.keyRelease(KeyEvent.VK_SHIFT);
%             rob.keyRelease(KeyEvent.VK_G);
% function [] = Setup_Solve_CarOCP(Position)
echo off
%This file was auto-exported via export_custom_fun.m with all necessary dependancies
mfile_name    = mfilename('fullpath');
[pathstr,~,~] = fileparts(mfile_name);
cd(pathstr);
% clearvars
addpath('./required_add_funs')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global state 
state = 1;

%%
% clear all
% close all
% clc
% clearvars -except xy_GI_ed xy_I1

import casadi.*

screen = 2;
fig_height = 100;
fig_width  = 50;
fig_pos    = 'middle-right';
% fig_pos = [Position(1) Position(2)];
% fig_height = Position(4);
% fig_width = Position(3);


fontsize = 18;
SetDefaultProperties("Fontsize",fontsize)
sim_CT = 1;
interp_method = 'pchip';

clearvars casadi_options ipopt_options
casadi_options = struct;
ipopt_options  = struct;
% ipopt
    ipopt_options.max_cpu_time      = 30;
    ipopt_options.max_iter          = 1000;
    ipopt_options.tol               = 10^-6;
    ipopt_options.dual_inf_tol      = 10^-6;
    ipopt_options.mumps_mem_percent = 16000; %8000;
%     ipopt_options.print_level       = 6;
%     ipopt_options.check_derivatives_for_naninf = 'yes';


%%
gray_bg_color = 0.9*[1 1 1];

trajectory_plots_on = 0;
animation_plots_on  = 1;
make_video_on       = 0;

close_IG_fig        = 1;
close_capture_fig   = 1;
close_video_fig     = 1;

plot_car_opt = 1;
load("optSolStruct.mat")


% video parameters
h_int_des = 10^-4;
fps = 10;
frametime = 1/fps;
h_int = (frametime/(ceil(frametime/h_int_des)*h_int_des)) * h_int_des; %rescale h s.t. N*h=1/fps
frameskips = round(frametime/h_int);
video_speed = 0.75;
window_size_coeff = 1/3;
xticks_h = 0.5;



%% problem specific I
IG_on = 1;
PathFollowing = 1;
avoid_objects = 1;
% load('matlab.mat')

max_abs_theta_2 = 0.4;

xy_0  = [0,-1];
xy_f  = [10,5];
psi_0 = 100*pi/180; %in rad, North
psi_f = 90*pi/180; %in rad, North

absb_delta = 40*pi/180; %in rad   %wheel angle
absb_omega = 90*pi/180; %in rad/s %wheel angle velocity
lb_v = -0.05;
ub_v = 8;    
add_xy = 5;
lb_xx  = -add_xy;
ub_xx  = 10+add_xy;
lb_y   = -add_xy;
ub_y   = 10+add_xy;
lb_a   = -10;
ub_a   = 5;


[rhs,n,m] = getCarDynamics;

%% problem specific II
%car
    [car,car_colors,indiv_part,circles,body_center] = create_car_object('blue');
    car_now = rot_transl_object(car,psi_0-pi/2,'rad',xy_0',body_center,indiv_part);

     fig = MakeDefaultFig(fig_width,fig_height,'Screen',screen,'FigPos',fig_pos,'FigTitle','Draw initial guess');
%     fig = app.UIAxes ;
    hold on
    patch_object(car_now,car_colors,1);
    my_rectangle(xy_f'-0.5,xy_f'+0.5);
    axis equal
    axis([lb_xx ub_xx lb_y ub_y])

% objects
    if avoid_objects==1
%         objects_circle_midpoint = [4;4];%[[4;4],[4;7]];
%         objects_cirlce_radius   = repmat(2.5,1,size(objects_circle_midpoint,2));
        if 0%isempty(objects_circle_midpoint) || isempty(objects_cirlce_radius)
            [objects_circle_midpoint,objects_cirlce_radius] = add_objects;
        else
            load('objects.mat')
        end
        for oo=1:size(objects_circle_midpoint,2)
    %         circle(objects_circle_midpoint(:,oo),2*objects_cirlce_radius(oo),'red');
            circle3(objects_circle_midpoint(:,oo),2*objects_cirlce_radius(oo),'FaceColor','red','FaceAlpha',0.8);
        end
    end

% object avoidance
    if avoid_objects==1
        safety_dist = 0.005;
        car_circle_midpoint = circles.midpoint;
        car_cirlce_radius   = circles.radius;
        g_oa = [];
        lb_g_oa = -Inf*ones(length(objects_cirlce_radius)*length(car_cirlce_radius),1);
        ub_g_oa =    0*ones(length(objects_cirlce_radius)*length(car_cirlce_radius),1);
        x_k = MX.sym('x_k', n);
        for oo=1:size(objects_circle_midpoint,2)
            for rr=1:length(car_cirlce_radius)
                car_circle_midpoint_now = x_k(1:2)+rot_z(car_circle_midpoint(:,rr),x_k(4)-pi/2,'rad');
                dist_vec = car_circle_midpoint_now - objects_circle_midpoint(:,oo);
                min_dist = car_cirlce_radius(rr)+objects_cirlce_radius(oo)+safety_dist;
                dist_fun = -(dist_vec(1)^2 + dist_vec(2)^2 - min_dist^2);
                g_oa = [g_oa; dist_fun]; 
            end
        end
        g_oa_fun = Function('g_oa_fun',{x_k},{g_oa});
    else
        g_oa_fun = [];
    end


%% SetOptsParams
% problem formulation                              
    OptsParams.vars_as_block   = 0; 
    OptsParams.constr_X_interm = 1;
    OptsParams.scale_x         = [];
    OptsParams.scale_u         = [];
    OptsParams.parametrize_OP  = 0;
    OptsParams.option_set      = 1;
    OptsParams.t_grid_normed_opt = 1;
    OptsParams.casadi_options = [];
% solver optionxs                       
    OptsParams.solver_choice  = 'ipopt'; %'ipopt','worhp'
    OptsParams.solver_options = ipopt_options; 
% integrator choice + cost functional
    integr_method_cell = {'ERK','IRK','IMP','ERK_adaptive','SPRK'};
    L_int_type_cell    = {'Mayer','GL','analytic'};
%_____________________________________________________________
    OptsParams.integr_method = integr_method_cell{1};
    OptsParams.L_int_type    = L_int_type_cell{1}; % analytic does not consider x_ki yet!
    OptsParams.deltaEMP_Mayer_on = 0;
% butcher tableau choice + GL choice
    des_int_order.ERK  = '4,1';
    des_int_order.IRK  = '3,1';
    des_int_order.ERK_adaptive  = '12,1';
    des_int_order.SPRK = '4,1';
        EMP_diff = 4*10^-4; %exponent has to be in {-1,-2,-3,...}
%     integr_method_params.eta_ERK2 = 0.5 + EMP_diff;
    integr_method_params.delta_MP = EMP_diff;
   %_____________________________________________________________
    OptsParams.des_int_order = des_int_order;
    OptsParams.integr_method_params = integr_method_params;       
    OptsParams.n_GL = 2;
% PATH FOLLOWING
    OptsParams.PathFollowing = PathFollowing;
% execute
    OptsParams      = removeEmptyFields(OptsParams);
    OptsParams_NVPs = namedargs2cell(OptsParams);
    OptsParams_WS   = A_SetOptsParams(OptsParams_NVPs{:});
    clearvars OptsParams_NVPs
 

%% Setup I
% choose T, h_des & M --> N                        
    QuantsOCP.T     = 15;      %             exact PH if QuantsOCP.w_J.t_f==0 and     (PH = prediction horizon)
                               % initial guess for PH if QuantsOCP.w_J.t_f>0
    QuantsOCP.h_des = 2*10^-1; % h<=h_des is as close as possible to h_des and satifies N*M*h = T (exactly)    (N = ceil(1/M * T/(h_des) and h=T/(N*M) <=> N*M*h=T)
    QuantsOCP.M     = 1;       % length of each single shooting interval   
    QuantsOCP.N     = [];%[];  % ~isempty(N)==1 --> N was chosen. This means that T, h_des & M are not involved in determining h.
% dynamics
    [rhs,n,m] = getCarDynamics;
    n = n;
    m = m;
% cost functionl
    R = diag([1 0.1]);
    Q = eye(n);
    QuantsOCP.L   = @(t,x,u,w_L) w_L*1/2*(u'*R*u + 0*x'*Q*x);
    QuantsOCP.Phi = [];
    QuantsOCP.y   = [];
% scale and weigh parts of cost functional   
    w_J.L   = 0;
    w_J.Phi = 0;
    w_J.t_f = 1;  
   %_____________________________________________________________
    QuantsOCP.w_J           = w_J; % if w_J.t_f>0 then t_f is part of J and is variable (in the below given bounds)!
    QuantsOCP.J_scale_coeff = [];  %[c_J_L, c_J_Phi, c_J_t_f]
% t_f                                              
    QuantsOCP.lb_t_f    = 10^-6; % lower bound on final time t_f
    QuantsOCP.ub_t_f    = 12;    % upper bound on final time t_f
    QuantsOCP.t_f_guess = 6;    % upper bound on final time t_f
% x    
    % IC and general bounds 
        QuantsOCP.x_0  = [xy_0 0 psi_0 0]';
        QuantsOCP.lb_x = [lb_xx  lb_y  lb_v -Inf -absb_delta]';
        QuantsOCP.ub_x = [ub_xx  ub_y  ub_v  Inf  absb_delta]';
    % constraint terminal state
        QuantsOCP.constr_x_f        = 1; 
        QuantsOCP.constr_states_x_f = [1,2];%,4];
        QuantsOCP.box_states        = 1:2;
        QuantsOCP.delta_box_x       = 0.5; %0 --> equality is forced
        QuantsOCP.x_f               = [xy_f 0 psi_f 0]';
    % guess
        QuantsOCP.x_guess_fun = [];
        QuantsOCP.x_guess_fun_further_params = [];
    % path constraints
        QuantsOCP.x_path_constr    = g_oa_fun; %Casadi Function of x_k    
        QuantsOCP.x_path_constr_lb = lb_g_oa;
        QuantsOCP.x_path_constr_ub = ub_g_oa;
% u                                                
    % general bounds
        QuantsOCP.lb_u = [lb_a  -absb_omega]';
        QuantsOCP.ub_u = [ub_a   absb_omega]';
    % constraint initial and terminal input
        QuantsOCP.constr_u_0  = 0;
        QuantsOCP.constr_u_f  = 0;
        QuantsOCP.delta_box_u = 0;
    % initial and terminal u
        QuantsOCP.u_0 = zeros(m,1);
        QuantsOCP.u_f = [];
    % guess
        QuantsOCP.u_guess_fun = [];
        QuantsOCP.u_guess_fun_further_params = [];
% du                                               
    QuantsOCP.lb_du_s = -[10; 300*pi/180];
    QuantsOCP.ub_du_s =  [10; 300*pi/180];
% y                                                
    QuantsOCP.lb_y = [];
    QuantsOCP.ub_y = [];
% scaling of state and input                                        
    QuantsOCP.T_x_scale = [];
    QuantsOCP.T_u_scale = [];

% PATH FOLLOWING
    % path-state dynamics
        QuantsOCP.rhs_path = @(t,theta,v) [50*theta(1) + v(1),-10*theta(2) + v(2)]; 
        QuantsOCP.n_p = 2;
        QuantsOCP.m_p = 2;
        QuantsOCP.xy_corridor = 'see below';
    % theta 
        %QuantsOCP.theta_0  = [xy_0 0 psi_0 0]';
        QuantsOCP.lb_theta = [0 -max_abs_theta_2]';
        QuantsOCP.ub_theta = [1  max_abs_theta_2]';
%         QuantsOCP.lb_theta = [-1 -max_abs_theta_2]';
%         QuantsOCP.ub_theta = [ 0  max_abs_theta_2]';
    % v                                              
        QuantsOCP.lb_v = [-Inf -Inf]';
        QuantsOCP.ub_v = [ Inf  Inf]';

% execute
    QuantsOCP      = removeEmptyFields(QuantsOCP);
    QuantsOCP_NVPs = namedargs2cell(QuantsOCP);
    SetupI_WS      = B_SetupI(rhs,n,m,OptsParams_WS.parametrize_OP,OptsParams_WS.PathFollowing,QuantsOCP_NVPs{:});
    clearvars QuantsOCP_NVPs

if exist('returnAfterSetupI')
    return
end

%% problem_initial guess
    if IG_on==1
        [xy_GI_ed,xy_I1,t_grid_ed_eval] = InitialGuessXY(xy_0,xy_f,'N',SetupI_WS.N,'interp_method',interp_method);
    end
    [xy_corridor,~,alpha_GI_ed] = xyCorridorParam(xy_I1,t_grid_ed_eval,interp_method);
    plotCorridor(xy_corridor,max_abs_theta_2)
    axes = gca;
    axes.Color = gray_bg_color;
    ca_children = get(gca,'children');
    set(gca,'children',ca_children([2:4,1,5:length(ca_children)]))
    SetupI_WS.xy_corridor = @(theta) xy_corridor(theta);
    drawnow

%     fig.WindowState = "minimized";
%     SetupI_WS.x_guess_norm_fun = @(t,x_0,x_f) [xy_GI_ed(t)'; x_0(3:n) + (x_f(3:n)-x_0(3:n))*t];
    SetupI_WS.x_guess_norm_fun = @(t,x_0,x_f) [xy_GI_ed(t)'; ub_v; -pi/2-alpha_GI_ed(t); 0];
%     SetupI_WS.x_guess_norm_fun = @(t,x_0,x_f) [xy_GI_ed(-1+t)'; ub_v; -pi/2-alpha_GI_ed(-1+t); 0];
echo on
fig(close);
save('Workspace.mat')
clear all
close all
clc
end

