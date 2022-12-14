function varargout = MakeDefaultFig(width,height,add)
    arguments
        % mandatory
            % nothing
        % optional 
            width double  = 30
            height double = 100
            add.FigTitle = []
            add.FigUnits  {mustBeMember(add.FigUnits, {'normalized','pixels'})} = 'normalized'
            add.DefaultAxesPosition double = [0.09, 0.1, 0.87, 0.84];
            add.Screen = 2
            add.FigPos {mustBeMember(add.FigPos, ...
                                     {'top-left',   'top-middle',   'top-right',      ...
                                      'middle-left','middle-middle','middle-right',   ...
                                      'bottom-left','bottom-middle','bottom-right'})} ...
                                        = 'middle-right';
            add.CallbackListenerFunOn = 1
            add.CallbackListenerFun   = @(src, evt) tightAxes
    end
    struct2CallerWS(add)

    fig_NVPs = {'units',FigUnits};
    if ~isempty(FigTitle)
        fig_NVPs = [fig_NVPs,'Name',FigTitle,'NumberTitle','off'];
    else   
        fig_NVPs = [fig_NVPs,'NumberTitle','on'];
    end
    if ~isempty(DefaultAxesPosition)
        fig_NVPs = [fig_NVPs,'DefaultAxesPosition',DefaultAxesPosition];
    end        

    fig = figure(fig_NVPs{:});
    
    if strcmp(FigUnits,'pixels')==1
        MonitorPositions = getScreenSize;
        if size(MonitorPositions,1)==1
            Screen = 1;
        end
        width  = width/MonitorPositions(Screen,3)*100;
        height = height/MonitorPositions(Screen,4)*100;
    end


%     cla
    setFigSize('Width',width,'Height',height,'FigPos',FigPos,'Screen',Screen,'FigHandle',fig)
    drawnow 

    
    if CallbackListenerFunOn==1 && ~isempty(CallbackListenerFun)
        hold on
        cur_axes = gca;
        listeneeer = setCallbackListenerFun(fig,cur_axes,CallbackListenerFun);
        CallbackListenerFun([],[])
    end
    
    
    if nargout==1
        varargout{1} = fig;
    end
end