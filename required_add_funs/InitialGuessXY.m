function [xy_GI_ed,xy_I1,t_grid_ed_eval] = InitialGuessXY(xy_0,xy_f,options)
    arguments
        % mandatory
            xy_0
            xy_f
        % optional 
            options.Fig = gcf;
            options.N = 10;
            options.N_fine_coeff = 20;
            options.interp_method = 'pchip';
            options.choose_method = 1;
    end
    struct2CallerWS(options)
    N_fine = N_fine_coeff*N;

    mind_dist = 0.2;

    theta_se = [0 1];
%     theta_se = [-1 0];

    %% setup
    global interm_points choose_points p_curr_guess
    choose_points = 1;
    switch choose_method
        case 1
            %%
            interm_points = [];

            set(Fig,'KeyPressFcn',@(src,event) KeyDown(src,event,xy_0,xy_f,N,N_fine,interp_method,theta_se))
            set(Fig,'WindowButtonDownFcn',  @(src,event) ImageClickCallback(src,event,gca,xy_0,xy_f,N,N_fine,interp_method,0,mind_dist,theta_se));
            set(Fig,'WindowButtonMotionFcn',@(src,event) ImageClickCallback(src,event,gca,xy_0,xy_f,N,N_fine,interp_method,1,mind_dist,theta_se));        
            p_curr_guess = plotCurrentGuess(interm_points,xy_0,xy_f,N,N_fine,interp_method,theta_se);
            drawnow
            while choose_points==1
                pause(eps)
            end
            set(Fig,'KeyPressFcn',[])
            set(Fig,'WindowButtonDownFcn',[])
            set(Fig,'WindowButtonMotionFcn',[])
    
        case 2
            %%
            interm_points = [];
            while choose_points==1
                cla        
                plotCurrentGuess(interm_points,xy_0,xy_f,N,N_fine);
            
                [x_interm,y_interm,button] = ginput(1);
                if button==1%==left mouse click
                    interm_points = [interm_points; [x_interm,y_interm]];
                elseif button==8%==backspace
                    if size(interm_points,1)>1
                        interm_points = interm_points(1:end-1,:);
                    elseif size(interm_points,1)==1
                        interm_points = [];
                    end
                elseif isempty(button)==1 %-> enter
                    choose_points = 0;
                end
            %     interm_points
            end
    end

    xy_GI = griddedInterpolant(theta_se(1):1/(2+size(interm_points,1)-1):theta_se(2),[xy_0; interm_points; xy_f],interp_method);
    xy_I1 =            interp1(theta_se(1):1/(2+size(interm_points,1)-1):theta_se(2),[xy_0; interm_points; xy_f],interp_method,'pp');

    t_grid_ed = compEquidistPoints(xy_GI,N,theta_se);
    t_grid_ed_eval = t_grid_ed(linspace(theta_se(1),theta_se(2),N_fine));
    xy_guesss_fine = xy_GI(t_grid_ed_eval);
    xy_GI_ed = griddedInterpolant(t_grid_ed_eval,xy_guesss_fine,interp_method);



end

%% auxiliary functions
function t_grid_ed = compEquidistPoints(xy_GI,N,theta_se)
    N_l     = 5*N;
    l_last  = Inf;
    while true
        xy = xy_GI(linspace(theta_se(1),theta_se(2),N_l));
        l_now = sum(vecnorm(diff(xy),2,2));
        delta_l = abs(l_now-l_last);
        if delta_l<1*10^-4
            break
        end       
        l_last = l_now;
        N_l = ceil(N_l*1.2);
    end
    % N_l
    
    l_cumsum = cumsum(vecnorm(diff(xy),2,2));
    help = l_cumsum/l_cumsum(end);
    t_grid_ed = griddedInterpolant([0; help],linspace(theta_se(1),theta_se(2),length(help)+1),'linear');
end

function varargout = plotCurrentGuess(interm_points,xy_0,xy_f,N,N_fine,interp_method,theta_se)
    xy_GI = griddedInterpolant(theta_se(1):1/(2+size(interm_points,1)-1):theta_se(2),[xy_0; interm_points; xy_f],interp_method);

    if 0
        t_grid_ed = compEquidistPoints(xy_GI,N);
    
        xy_guesss      = xy_GI(t_grid_ed(linspace(theta_se(1),theta_se(2),N)));
        xy_guesss_fine = xy_GI(t_grid_ed(linspace(theta_se(1),theta_se(2),N_fine)));
    
        p_curr_guess(1) = plot(xy_guesss_fine(:,1),xy_guesss_fine(:,2),'b-','LineWidth',1.5);
        try
            p_curr_guess(2) = plot(xy_guesss(:,1),xy_guesss(:,2),        'ko','LineWidth',2,'MarkerSize',5);
            p_curr_guess(3) = plot(interm_points(:,1),interm_points(:,2),'r+','LineWidth',2,'MarkerSize',10); %fails for interm_points==[]
        end
    else
        xy_guesss_fine = xy_GI(linspace(theta_se(1),theta_se(2),N_fine));
        p_curr_guess(1) = plot(xy_guesss_fine(:,1),xy_guesss_fine(:,2),'b-','LineWidth',1.5);
        try
            p_curr_guess(2) = plot(interm_points(:,1),interm_points(:,2),'r+','LineWidth',2,'MarkerSize',10); %fails for interm_points==[]
        end
    end


    varargout{1} = p_curr_guess;
end

%%
function KeyDown(src,event,xy_0,xy_f,N,N_fine,interp_method,theta_se)
    global interm_points choose_points p_curr_guess
    switch event.Key
        case 'backspace'
            if size(interm_points,1)>1
                interm_points = interm_points(1:end-1,:);
            elseif size(interm_points,1)==1
                interm_points = [];
            end
        case 'return'
            if ~isempty(interm_points)            
                choose_points = 0;
            else
                warning('choose at least one point')
            end
    end
    delete(p_curr_guess)
    p_curr_guess = plotCurrentGuess(interm_points,xy_0,xy_f,N,N_fine,interp_method,theta_se);
end

function ImageClickCallback(src,event,hAxes,xy_0,xy_f,N,N_fine,interp_method,hover,mind_dist,theta_se)
    global interm_points p_curr_guess
    coordinates = get(hAxes,'CurrentPoint');
    coordinates = coordinates(1,1:2);

    if isempty(interm_points)==1 || norm(interm_points(end,:)-coordinates)>=mind_dist
        interm_points_plot = [interm_points; coordinates];
        if hover==0
            interm_points = [interm_points; coordinates];
        end
    else
        interm_points_plot = interm_points;
    end

    delete(p_curr_guess)
    p_curr_guess = plotCurrentGuess(interm_points_plot,xy_0,xy_f,N,N_fine,interp_method,theta_se);
end