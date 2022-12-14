function listeneeer = setCallbackListenerFun(cur_fig,cur_axes,listener_fun)
    arguments
        % mandatory
            % nothing
        % optional 
            cur_fig  = gcf;
            cur_axes = gca;
            listener_fun = @(src, evt) tightAxes
    end

    %% change of strings
    string_props = {'XLabel','YLabel','ZLabel','Title'};
    string_props_array = [];
    for ii=1:length(string_props)
        string_props_array = [string_props_array, cur_axes.(string_props{ii})];
    end
    listeneeer{1} = addlistener(string_props_array, 'String', 'PostSet', listener_fun);
    
    %% axes
    if 1
        fns_axes = fieldnames(cur_axes);
        not_inds = [];
        excl_props = {'Position','CurrentPoint','TightInset','Children','Layout','Type','Legend','HandleVisibility'};
        for ii=1:length(fns_axes)
            if any(strcmp(fns_axes{ii},excl_props))==1
                not_inds = [not_inds, ii];
            end
        end
        prop_axes = fns_axes(setdiff(1:end,not_inds));
        listeneeer{2} = addlistener(cur_axes, prop_axes, 'PostSet', listener_fun);
    end
    
    %% 
%     listeneeer{3} = event.listener(cur_fig, 'LocationChanged', listener_fun);
    cur_fig.SizeChangedFcn = listener_fun;      


%     listener(cur_fig,'CurrentAxes','PostSet', listener_fun);
    
    %% figure
    if 0
    fns_fig = fieldnames(cur_fig);
    not_inds = [];
    for ii=1:length(fns_fig)
        if any(strcmp(fns_fig{ii},{'Position','OuterPosition','InnerPosition','Children','Scrollable','Type','WindowState'}))==1
            not_inds = [not_inds, ii];
        end
    end
    prop_fig = fns_fig(setdiff(1:end,not_inds));
    addlistener(cur_fig, prop_fig, 'PostSet', listener_fun)
    end

end