function [] = second()
warning('off','all');
echo off
load("Workspace.mat");
addpath('.\required_add_funs');

disp('initialisiere Optimizer');
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
% clc
% return
% x_traj_CT_GI  = griddedInterpolant(tgrid_x_CT,x_traj_CT','spline');
% optSol = var2struct(x_traj_CT_GI,t_f_opt);

window_size = window_size_coeff*t_f_opt;

car_colors_opt = car_colors;
car_colors_opt(1:2) = {[0 1 0],[0 1 0]};
if animation_plots_on==1    
    x_0_guess = full(SolveOCP_WS.sol.x);
    if avoid_objects==1
        x_lim_help = [x_traj_CT(1,:),x_0_guess(1,:),objects_circle_midpoint(1,:)];
        y_lim_help = [x_traj_CT(2,:),x_0_guess(2,:),objects_circle_midpoint(2,:)];
    else
        x_lim_help = [x_traj_CT(1,:),x_0_guess(1,:)];
        y_lim_help = [x_traj_CT(2,:),x_0_guess(2,:)];
    end

    fig_capture_frames = MakeDefaultFig(fig_width,fig_height,'CallbackListenerFunOn',0, ...
        'FigPos',fig_pos,'Screen',screen,'FigTitle','Capture frames for video');
    jFrame = get(handle(fig_capture_frames), 'JavaFrame');
    jFrame.setMinimized(1);
    %     set(gcf,'WindowState','minimized')
    p = panel();
    p.pack({1}, {2/3 []});
    p(1,2).pack(2,1);
%         height_p12 = 1/2+0.05;
%         p(1,2).pack({height_p12/2 height_p12/2}, {1});
    p.select('all');
    p.margin = 15;
    p.margintop = 8;
    p.marginbottom = 20;
    p.marginleft = 18;
    p(1,2,1).marginbottom = 5;
    p(1,1).marginright = 27;
    p(1,2,1,1).select();
%         ylabel('input $u$')
        ylabel('Eingangsgr{\"o}{\ss}en')
        set(gca, 'xtick', []);
    p(1,2,2,1).select();
%         xlabel('time $t$ in seconds')
%         ylabel('parts of the state $x$')
        xlabel('Zeit $t$ in Sekunden')
        ylabel('Fahrzeugzust{\"a}nde')
        xticks(0:xticks_h:t_f_opt)
    p(1,1).select();
%         xlabel('$x$-position')
%         ylabel('$y$-position')
        xlabel('Position in $x$-Richtung')
        ylabel('Position in $y$-Richtung')
        hold on
        my_rectangle(SetupII_WS.lb_x_f([1,2]),SetupII_WS.ub_x_f([1,2]));
        add_space = 1.5;
        axis equal
        set(gca,'Color',gray_bg_color)
        plotCorridor(xy_corridor,max_abs_theta_2)
        plot(x_traj_CT(1,:),x_traj_CT(2,:),'k-','Linewidth',2)
%         plot(x_traj(1,:),x_traj(2,:),'mo','Linewidth',2,'LineStyle','none')
        scatter(x_traj(1,:),x_traj(2,:),'MarkerFaceColor','none', ...
            'MarkerEdgeColor',0.5*[1 1 1],'MarkerEdgeAlpha',0.8)
        plot(x_0_guess(1,:),x_0_guess(2,:),'k--')  
        if avoid_objects==1
            for oo=1:length(objects_cirlce_radius)
    %             circle(objects_circle_midpoint(:,oo),2*objects_cirlce_radius(oo),'red');
    %             all_circles(oo) = 
                circle3(objects_circle_midpoint(:,oo),2*objects_cirlce_radius(oo),'FaceColor','red','FaceAlpha',1);
            end
        end
        xlim([ min(x_lim_help)-add_space,...
               max(x_lim_help)+add_space ])
        ylim([ min(y_lim_help)-add_space,...
               max(y_lim_help)+add_space ])

    drawnow
    tgrid_x_CT_int = 0:h_int:t_f_opt;
    x_traj_CT_int  = interp1(tgrid_x_CT,x_traj_CT',tgrid_x_CT_int,'spline')';
    p(1,2,1,1).select();
        hold on
        plot(tgrid_x_CT_int,u_CT(tgrid_x_CT_int),'LineWidth',2)
%         plot(tgrid_x_CT_int(1:end-1),diff(u_CT(tgrid_x_CT_int)),'LineWidth',2)
        set_limits_perc(u_CT(tgrid_x_CT_int),[8,8]);
%         leg_121 = legend('$a$','$\omega$','Location','east');
        leg_121 = legend('Beschleunigung','Lenkwinkelgeschw.','Location','southeast');
        setLegColorAlpha(leg_121,[1 1 1],0.7)
    p(1,2,2,1).select();
        hold on
        plot(tgrid_x_CT_int,x_traj_CT_int(3:end,:),'LineWidth',2)
        set_limits_perc(x_traj_CT_int(3:end,:),[8,8]);
        leg_122 = legend('Geschwindigkeit','Fahrzeugwinkel (rad)','Radwinkel (rad)','Location','east');
        setLegColorAlpha(leg_122,[1 1 1],0.7)
%         [car,car_colors,indiv_part,circles,body_center] = create_car_object;
    clearvars alpha frames frames_opt

    t_plot    = 0;
    t_frames  = 0;

    frame_numbers = [1:frameskips:length(tgrid_x_CT_int), length(tgrid_x_CT_int)];

    t_anim_0 = datevec(datetime('now'));
    f = waitbar(0,'Rendering Video');
%     f.Position= [900.2500 624 547.5000 90.5000];
for jj=1:length(frame_numbers)
        waitbar(jj/length(frame_numbers),f,'Rendering Video');
        ii = frame_numbers(jj);
        tic
        p(1,1).select();
            x_i = x_traj_CT_int(:,ii);
            indiv_part.angle = x_i(5)*ones(1,2); 
            car_now = rot_transl_object(car,x_i(4)-pi/2,'rad',x_i([1,2]),body_center,indiv_part);
            all_car_patches = patch_object(car_now,car_colors,0.65);
            dot = plot(x_i(1),x_i(2),'yo','Linewidth',3);    
            for rr=1:length(car_cirlce_radius)
                car_circle_midpoint_now = x_i(1:2)+rot_z(car_circle_midpoint(:,rr),x_i(4)-pi/2,'rad');
                car_circles_now(rr) = circle2(car_circle_midpoint_now,2*car_cirlce_radius(rr),'k');                    
            end

            if ii~=frame_numbers(end)
                title(sprintf('t = %0.2f sek',tgrid_x_CT_int(ii)),'interpreter','none')
            else
                title(sprintf('t = %0.5f sek',tgrid_x_CT_int(ii)),'interpreter','none')
%                 title(sprintf('$t = %0.5f \\ \\mathrm{sek}$',tgrid_x_CT_int(ii)),'interpreter','tex')
%                 title(sprintf('$t = %0.5f \\ \\mathrm{sek}$',tgrid_x_CT_int(ii)),'interpreter','tex')
            end
        
        p(1,2,1,1).select();
            xlim([max(0, tgrid_x_CT_int(ii)-window_size/2), min(tgrid_x_CT_int(ii)+window_size/2, tgrid_x_CT_int(end))])
            dots_1 = plot(tgrid_x_CT_int(ii),u_CT(tgrid_x_CT_int(ii)),'k-o','LineWidth',2);
        p(1,2,2,1).select();
            xlim([max(0, tgrid_x_CT_int(ii)-window_size/2), min(tgrid_x_CT_int(ii)+window_size/2, tgrid_x_CT_int(end))])
            dots_2 = plot(tgrid_x_CT_int(ii),x_traj_CT_int(3:end,ii),'k-o','LineWidth',2);
        drawnow
        t_plot = t_plot+toc;

        tic
        if 0%jj>1 && plot_car_opt==1 && t_now>optSol.t_f_opt
            for mm=1:length(all_car_patches_opt)
                all_car_patches_opt(mm).Visible = 'off';
            end
        end
        cur_frame = getframe(gcf);
        if 0%jj>1 && plot_car_opt==1 && t_now>optSol.t_f_opt
            for mm=1:length(all_car_patches_opt)
                all_car_patches_opt(mm).Visible = 'on';
            end
        end
        frames(jj) = cur_frame;
        t_frames = t_frames + toc;

        t_now = tgrid_x_CT_int(ii);
        if plot_car_opt==1 
            if 1%t_now<=optSol.t_f_opt
                p(1,1).select();
                if t_now<=optSol.t_f_opt
                    t_eval = t_now;
                else
                    t_eval = optSol.t_f_opt;
                end
                x_i = optSol.x_traj_CT_GI(t_eval)';
                indiv_part.angle = x_i(5)*ones(1,2); 
                car_now = rot_transl_object(car,x_i(4)-pi/2,'rad',x_i([1,2]),body_center,indiv_part);
                all_car_patches_opt = patch_object(car_now,car_colors_opt,0.2);
            end
            cur_frame = getframe(gcf);
            frames_opt(jj) = cur_frame;
        end

        if jj<length(frame_numbers)
            delete(all_car_patches);
            delete(dot)
            delete(dots_1)
            delete(dots_2)
            delete(car_circles_now)
            t_next = tgrid_x_CT_int(frame_numbers(jj+1));
            if 1%plot_car_opt==1 &&  t_next<=optSol.t_f_opt
                delete(all_car_patches_opt)
            end
        end
    end
    t_anim_f = datevec(datetime('now'));
    t_anim   = etime(t_anim_f,t_anim_0);
    [t_anim_h,t_anim_m,t_anim_s] = hms(seconds(t_anim));
%     fprintf('Animation took %0.1fsec\n',t_anim_s)

    PrintFig('final_frame','FileFormat','pdf')
    if close_capture_fig==1
        close(fig_capture_frames)
    end

     name_person = 'Markus';
     video_title = [sprintf('%0.5fs',t_f_opt),'_',name_person,'_at_',datestr(now,'hhMM')];
    if PathFollowing==1
        video_title = [video_title,'_PF'];
    end

%     fig_video = MakeDefaultFig(fig_width,fig_height,'CallbackListenerFunOn',0,'Screen',screen,'FigTitle','Video','FigPos',fig_pos);
%     set(gcf,'DefaultAxesPosition',[0 0 1 1])
%     movie(frames,1,round(video_speed*fps))
%     movie(frames_opt,1,round(video_speed*fps))
%     if close_video_fig==1
%         close(fig_video)
%     end
    if make_video_on==1
        framerate = ceil(video_speed*fps);
        make_video([frames,repmat(frames(end),1,framerate+1)],framerate,video_title);
    end
end

close(f);
save('Workspace.mat','fig_height','fig_width','fig_pos','frames_opt','video_speed','fps','screen','t_f_opt');

clear all
close all


end

function dx = eval_rhs(t,x,u,rhs)
    dx = rhs(t,x,u(t)');
end
