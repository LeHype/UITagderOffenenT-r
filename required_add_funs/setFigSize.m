function setFigSize(options)
    arguments
        % mandatory
            % nothing
        % optional 
            options.FigHandle                  = gcf;
            options.Width         (1,1) double = 50;
            options.Height        (1,1) double = 100;
            options.TaskbarHeight (1,1) double = 42;
            options.FigPos        {mustBeMember(options.FigPos, ...
                                        {'top-left',   'top-middle',   'top-right',      ...
                                         'middle-left','middle-middle','middle-right',   ...
                                         'bottom-left','bottom-middle','bottom-right'})} ...
                                               = 'middle-right';
            options.Screen {mustBeNumeric}     = 1;         
            options.PositionName {mustBeMember(options.PositionName,{'InnerPosition','OuterPosition'})} ...
                                               = 'OuterPosition';
    end
    struct2CallerWS(options)

    set(0,'units','pixels');  
    pixel_all_displays = get(0,'MonitorPositions'); % does not work if monitor layout changes within one Matlab session
%     [~,displays_sort_ind]  = sort(pixel_all_displays(:,1));
    displays_sort_ind = 1:size(pixel_all_displays,1);
    pixel_all_displays = pixel_all_displays(displays_sort_ind,:);
    n_Screens = size(pixel_all_displays,1);
    if n_Screens==1 && Screen~=1 %only one screen is available --> "overwrites" Screen 
        Screen = 1;
    elseif Screen>n_Screens
        Screen = 1;
    end
    pixel = pixel_all_displays(Screen,:);
    pix_w = pixel(3);  
    pix_h = pixel(4) - TaskbarHeight;
    l_monitor_corr = pixel(1);
    b_monitor_corr = pixel(2);
    b_f_help = pix_h - (Height/100)*pix_h; 
    l_f_help = pix_w - (Width/100)*pix_w;
%     if nargin>=5
        vert = FigPos(1:find(FigPos=='-')-1);
        hor  = FigPos(find(FigPos=='-')+1:length(FigPos));
        if strcmp(vert,'bottom')
            b_f = b_monitor_corr + 0*b_f_help;% + taskbar_h;
        elseif strcmp(vert,'middle')
            b_f = b_monitor_corr + 0.5*b_f_help;% + taskbar_h;
        elseif strcmp(vert,'top')
            b_f = b_monitor_corr + 1*b_f_help;% + taskbar_h;  
        else
            error('Choose proper (vertical) position')
        end
        if strcmp(hor,'left')
            l_f = l_monitor_corr + 0*l_f_help;
        elseif strcmp(hor,'middle')
            l_f = l_monitor_corr + 0.5*l_f_help;
        elseif strcmp(hor,'right')
            l_f = l_monitor_corr + 1*l_f_help;  
        else
            error('Choose proper (horizontal) position')
        end
%     else
%         l_f = l_f_help;
%         b_f = b_f_help ;
%     end
    w_f = (Width/100)*pix_w;  
    h_f = (Height/100)*pix_h;
    set(FigHandle,'units','pixels',PositionName,[l_f b_f+TaskbarHeight w_f h_f])
    set(FigHandle,'units','normalized')
end
