function [] = second(name)
warning('off','all');
echo off
load("Workspace.mat");
addpath('.\required_add_funs');

disp('Initialisiere Solver: Optimierung wird gestartet');
%SECOND Summary of this function goes here
%   Detailed explanation goes here
%% SetupII, SingleShootIntEval & BuildOCP 
SetupII_WS            = C_SetupII(OptsParams_WS,SetupI_WS);
SingleShootIntEval_WS = D_SingleShootIntEval(OptsParams_WS,SetupI_WS,SetupII_WS);
BuildOCP_WS           = E_BuildOCP(SingleShootIntEval_WS,OptsParams_WS,SetupI_WS,SetupII_WS);

% StructureOCPHessian


%% SolveAndPostProcessOCP
% try 
%     BuildOCP_WS.w_0 = full(SolveOCP_WS.sol.x);
% catch
%     warning('not successful')
% end
SolveOCP_WS = F_SolveOCP(BuildOCP_WS,OptsParams_WS,n,m);%,'casadi_options',casadi_options);
PostProcessOCP_WS = G_PostProcessOCP2(SolveOCP_WS.sol,SingleShootIntEval_WS.X_mat,BuildOCP_WS,OptsParams_WS,n,m,SetupI_WS.M,SetupI_WS.N);%, ...
    %xy_corridor);


%% CT simulation ode45
% struct2CallerWS(PostProcessOCP_WS)
struct2CallerWS(PostProcessOCP_WS,{'tgrid_u','tgrid_x','u_traj','x_traj','t_f_opt','plot_u'})
t_f_opt;

clearvars t_ode45 x_ode45
if sim_CT==1
    x_0 = QuantsOCP.x_0;
    u_CT_interp_method = 'previous';
    u_CT    = griddedInterpolant(tgrid_u,u_traj',u_CT_interp_method);
%     u_CT_pp =            interp1(tgrid_u,u_traj',u_CT_interp_method,'pp');
    [tgrid_x_CT,x_traj_CT] = ode45(@(t,x) eval_rhs(t,x,u_CT,rhs),[0 t_f_opt], x_0, ...
                                              odeset('RelTol',10^-13,'AbsTol',10^-16));
    x_traj_CT = x_traj_CT';
end

if SolveOCP_WS.solver.stats.success;
disp('------------------------------')
    disp('Optimieren erfolgreich');
disp('------------------------------')
else 
    warning('on','all');
    warning('Optimieren Fehlgeschlagen! ')
    warning('off','all');
% disp('Optimieren Fehlgeschlagen! ')
disp('Dies kann merere Gründe haben.')
disp('  1. Der gezeichnete Pfad physikalisch ist nicht fahrbar.')
disp('     Habt ihr eventuell zu steile Kurven eingebaut? Und denkt ihr')
disp('     dass ein echtes Auto diesen Pfand fahren könnte?')
disp('  2. Der Optimizer ist in ein sog. lokales Minimum gelaufen.')
disp('     Sollte dies der Fall sein kann durch Anpassung des  ')
disp('     Algorythmus meistens die Lösung noch gefunden werden.')
disp('  3. Die maximale Fahrtzeit von 12 Sekunden wurde überschritten.')
disp('     Falls ihr einen besonders langen Pfad ohne zu steile Kuven ')
disp('     ausgewählt habt ist dies wahrscheinlich der Grund.')
disp('  4. Falls ihr euch immernoch nicht sicher seit was fehlgeschlagen ist')
disp('     fragt gerne eine Person vom Lehrstuhl.')
end
disp('------------------------------')
disp('Video der Fahrt wird nun gerendert')
disp('------------------------------')

%% trajectory plots
if trajectory_plots_on==1
    left_color  = [0 0 1];  %red
    right_color = [1 0 0];  %blue
    markertype = '.';
    
    plotted_states = 1:n;
    
    MakeDefaultFig(50,100,'Screen',screen,'CallbackListenerFunOn',0,'FigPos','middle-right', ...
        'FigTitle','Trajectories');
    p = panel();
    p.pack({1}, {1/2 1/2});
    % p.pack({1}, {1});
    n_plots_vert = 2;
    p(1,1).pack(num2cell(1/n_plots_vert*ones(1,n_plots_vert)), {1});
    p(1,1).marginright = 30;
    p.select('all');
    p.fontsize = fontsize;
    leg_text = {'$x$','$y$','$v$','$\psi$','$\delta$'};
    
    % p(1,1,1,1) Discrete Time
        p(1,1,1,1).select();
        hold on
        pl = stairs(tgrid_x,x_traj(plotted_states,:)','Marker',markertype);
    %     pl = plot(tgrid_x(1:end-1),diff(x_traj(plotted_states,:),1,2)','Marker',markertype);
        xlim([0 tgrid_u(end)])
        legend(leg_text{plotted_states},'Location','NorthWest')
    
    % p(1,1,2,1) Continuous Time
    try
        p(1,1,2,1).select();
            plot(tgrid_x_CT,x_traj_CT(plotted_states,:))
            xlim([tgrid_x_CT(1) tgrid_x_CT(end)])
            legend(leg_text{plotted_states},'Location','NorthWest')
            yline(QuantsOCP.x_f(1))
            yline(1/10*QuantsOCP.x_f(4))
        x_traj(SetupI_WS.constr_states_x_f,end)-x_traj_CT(SetupI_WS.constr_states_x_f,end)
    end
    
    % p(1,2) Input Trajectory
        p(1,2).select();
        hold on
        plot_u(tgrid_u,u_traj','Marker',markertype)
        xlim([0 tgrid_u(end)])
        legend('$u_r$','$u_{\theta}$','Location','East')
end

%%
%%
if close_IG_fig==1
    close(fig)
else
    figure(fig)
    plot(x_traj_CT(1,:),x_traj_CT(2,:),'k-','Linewidth',2)
    plot(x_traj(1,:),x_traj(2,:),'g--','Linewidth',2,'Marker','.','MarkerSize',20)
    plot(x_traj(1,:),x_traj(2,:),'g','Linewidth',2,'Marker','o','MarkerSize',5,'Linestyle','none')
end

%% animation

% x_traj_CT_GI  = griddedInterpolant(tgrid_x_CT,x_traj_CT','spline');
% optSol = var2struct(x_traj_CT_GI,t_f_opt);


if animation_plots_on==1    
    carHeight = 0.5;
    [dims,coords] = create3DCarObjectTdoT(1.8,90,'BodyHeight',carHeight,'BodyWidth',1,'FrontLength',0.5,'TireRadius',0.2,'TireWidth',0.1);
    coords.YellowCylinder = CoordsCylinder(0.1,carHeight,11,'Translation',[0 0 carHeight/2]');
    
    tgrid_x_int = 0:h_int:tgrid_x(end);
    
    x_fun_help = griddedInterpolant(tgrid_x,x_traj','spline');
    x_fun   = @(t) x_fun_help(t)';
    phi_fun = computeTireAngleFun(x_fun,tgrid_x_int);
    
    x_fun_opt   = @(t) optSol.x_traj_CT_GI(min(t,optSol.t_f_opt))';
    phi_fun_opt = computeTireAngleFun(x_fun_opt,0:h_int:optSol.t_f_opt);

    
    obstacles.radius   = objects_cirlce_radius;
    obstacles.midpoint = objects_circle_midpoint;


    figureData.screen = 2;
    figureData.height = 100;
    figureData.width = 90.5;
    figureData.position = 'middle-right';
    figureData.backgrColor = [0.9000 0.9000 0.9000];
    figureData.fontsize = 18;
    figureData.xticks_h = 0.5000;


    VideoSpeed = 0.75;
    VideoTitle = [num2str(t_f_opt),'_s_Fahrtzeit_',name];
    [vidObj,t_grid_video,VideoFPS_actual] = prepareVideoObj(tgrid_x,'VideoFPS',12,'VideoSpeed',VideoSpeed, ...
        'VideoFormat','MPEG-4','VideoTitle',VideoTitle);
    frames_opt = animate3DCarTdoT(dims,coords,obstacles,figureData,SetupII_WS,SolveOCP_WS, ...
        t_grid_video,tgrid_x_int,x_fun,x_fun_opt,u_CT,phi_fun,phi_fun_opt,xy_corridor,max_abs_theta_2,'CaptureFrames',1);

    PrintFig('final_frame','FileFormat','pdf')
    if close_capture_fig==1
        close(gcf)
    end

end

save('Workspace.mat','figureData','frames_opt','vidObj','VideoFPS_actual','VideoSpeed');

clear all
close all


end

function dx = eval_rhs(t,x,u,rhs)
    dx = rhs(t,x,u(t)');
end
