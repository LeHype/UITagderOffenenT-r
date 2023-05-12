function FRAMES = animate3DCarTdoT(dims,coords,obstacles,figureData,SetupII_WS,SolveOCP_WS, ...
    t_grid,t_grid_int,x_fun,x_fun_opt,u_fun,phi_fun,phi_fun_opt,xy_corridor,max_abs_theta_2,options)
arguments
    dims
    coords
    obstacles
    figureData
    SetupII_WS
    SolveOCP_WS
    t_grid
    t_grid_int
    x_fun
    x_fun_opt
    u_fun
    phi_fun
    phi_fun_opt
    xy_corridor
    max_abs_theta_2
    options.LastStateOnly {mustBeMember(options.LastStateOnly,[0,1])} = 0
    options.CaptureFrames {mustBeMember(options.CaptureFrames,[0,1])} = 0
    options.Directory     = []
    options.VideoTitle    = []
    options.FigureWidth   = 100
    options.FigureHeight  = 100
    options.FigureUnits  {mustBeMember(options.FigureUnits,{'normalized','pixels'})} = 'normalized'
    options.FigurePos     = 'middle-right'
    options.Screen        = 1
end
struct2CallerWS(options)

obstacleHeight = 1;
obstacles.midpoint=[obstacles.midpoint;zeros(1,size(obstacles.midpoint,2))];

window_size = 0.333*t_grid(end);


%% setup figure and subplots
x_0_guess = full(SolveOCP_WS.sol.x);
x_traj = x_fun(t_grid);
x_traj_int = x_fun(t_grid_int);

x_lim_help = [x_traj(1,:),x_0_guess(1,:),obstacles.midpoint(1,:)];
y_lim_help = [x_traj(2,:),x_0_guess(2,:),obstacles.midpoint(2,:)];

fig_capture_frames = MakeDefaultFig(figureData.width,figureData.height,'CallbackListenerFunOn',0, ...
    'FigPos',figureData.position,'Screen',figureData.screen,'FigTitle',VideoTitle);
jFrame = get(handle(fig_capture_frames), 'JavaFrame');
jFrame.setMinimized(1);

p = panel();
p.pack({1}, {2/3 []});
p(1,2).pack(2,1);
p.select('all');
p.margin = 15;
p.margintop = 8;
p.marginbottom = 20;
p.marginleft = 3;
p(1,2,1).marginbottom = 5;
p(1,1).marginright = 45;
p(1,2,1,1).select();
    ylabel('Eingangsgr{\"o}{\ss}en')
    set(gca, 'xtick', []);
p(1,2,2,1).select();
    xlabel('Zeit $t$ in Sekunden')
    ylabel('Fahrzeugzust{\"a}nde')
    xticks(0:figureData.xticks_h:t_grid(end))
p(1,1).select();
    xlabel('Position in $x$-Richtung')
    ylabel('Position in $y$-Richtung')
    view([7 32.3])
    % view([0 90])
    hold on
    my_rectangle(SetupII_WS.lb_x_f([1,2]),SetupII_WS.ub_x_f([1,2]));
    axis equal
    set(gca,'Color',figureData.backgrColor,'ZColor','none');
    plotCorridor(xy_corridor,max_abs_theta_2)
    plot(x_traj_int(1,:),x_traj_int(2,:),'k-','Linewidth',2)
    scatter(x_traj(1,:),x_traj(2,:),'MarkerFaceColor','none', ...
        'MarkerEdgeColor',0.5*[1 1 1],'MarkerEdgeAlpha',0.8)
    plot(x_0_guess(1,:),x_0_guess(2,:),'k--')  
    coords_Cylinder = CoordsCylinder(obstacles.radius(1),obstacleHeight,20);
    if 1%options.collisionAvoidance==1
        for oo=1:length(obstacles.radius)
            PatchCylinder(coords_Cylinder, 'Translation', obstacles.midpoint(:,oo),'Rotation',eye(3),'RotCenterShift',[0;0;-0.5])
        end
    end
    add_space = 1.5;
    xlim([ min(x_lim_help)-add_space,...
           max(x_lim_help)+add_space ])
    ylim([ min(y_lim_help)-add_space,...
           max(y_lim_help)+add_space ])
    % xlim([-1 3])
    % ylim([-1 3])
    zlim([0 obstacleHeight+0.5])
    grid on
    animAxes = gca;
% p(2).select();
%      copyaxes(animAxes,gca,true)


%% plot trajectories
p(1,2,1,1).select();
    hold on
    p1 = plot(t_grid_int,u_fun(t_grid_int),'LineWidth',2);
    set_limits_perc(u_fun(t_grid_int),[8,8]);
    leg_121 = legend(p1,'Beschleunigung','Lenkwinkelgeschw.','Location','southeast','Fontsize',figureData.fontsize-5,'AutoUpdate','off');
    setLegColorAlpha(leg_121,'LegColor',[1 1 1],'LegAlpha',0.7)
p(1,2,2,1).select();
    hold on
    p2 = plot(t_grid_int,x_traj_int(3:end,:),'LineWidth',2);
    set_limits_perc(x_traj_int(3:end,:),[8,8]);
    leg_122 = legend(p2,'Geschwindigkeit','Fahrzeugwinkel (rad)','Radwinkel (rad)','Location','east','Fontsize',figureData.fontsize-5,'AutoUpdate','off');
    setLegColorAlpha(leg_122,'LegColor',[1 1 1],'LegAlpha',0.7)
clearvars alpha frames frames_opt


%% for loop over all desired frames
uds_blue = 1/250*[0,72,119];
uds_red  = 1/250*[200,34,84];
green    = [0 1 0]; 

f = waitbar(0,'Rendering Video');
for ii=1:length(t_grid)
    waitbar(ii/length(t_grid),f,'Rendering Video');
    t_i = t_grid(ii);
    tic
    p(1,1).select();
        x_i = x_fun(t_i);
        phi_i = phi_fun(t_i);
        allCarPatches = plot3DCarTdoT(coords,dims,x_i,phi_i,uds_blue,uds_blue);        
        if ii~=length(t_grid)
            title(sprintf('t = %0.2f sek',t_i),'interpreter','none')
        else
            title(sprintf('t = %0.5f sek',t_i),'interpreter','none')
        end        
    p(1,2,1,1).select();
        xlim([max(0, t_i-window_size/2), min(t_i+window_size/2, t_grid(end))])
        dots_1 = plot(t_i,u_fun(t_i),'k-o','LineWidth',2);
    p(1,2,2,1).select();
        xlim([max(0, t_i-window_size/2), min(t_i+window_size/2, t_grid(end))])
        dots_2 = plot(t_i,x_i(3:end),'k-o','LineWidth',2);

    if 1%options.plotOptimalCar==1 
        p(1,1).select();
        x_i = x_fun_opt(t_i);  
        phi_i = phi_fun_opt(t_i);
        allCarPatchesOpt = plot3DCarTdoT(coords,dims,x_i,phi_i,green,green,0);
    end
    %getting the frames:
    if CaptureFrames==1
        FRAMES(ii) = getframe(gcf);
    else
        FRAMES = [];
    end
    drawnow
    if ii<length(t_grid)
        delete(allCarPatches);
        delete(dots_1)
        delete(dots_2)
        delete(allCarPatchesOpt)
    end

%     fprintf('Animation took %0.1fsec\n',t_anim_s)
end
close(f);


end