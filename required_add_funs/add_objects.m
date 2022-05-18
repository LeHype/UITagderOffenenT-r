function [objects_circle_midpoint,objects_cirlce_radius] = add_objects(options)
    arguments
        % mandatory  
            % nothing
        % optional 
            options.d = 1
            options.Fig = gcf
            options.Axes = gca
    end

    d = options.d;
    Fig = options.Fig;
    Axes = options.Axes;

    last_click_loc = [Inf,Inf];
    cur_click_loc  = [Inf,Inf];
    stop_exec = 0;
    not_done  = 1;
    
    set(Fig,'WindowButtonDownFcn',@MouseClickCallback);  
    set(Fig,'WindowKeyPressFcn',  @KeyPressCallback);
    set(Fig,'CurrentAxes', Axes)
    
    objects_circle_midpoint = [];
    while stop_exec==0 && not_done==1 
        pause(10^-2)
        if all(last_click_loc~=cur_click_loc)
            set(Fig, 'CurrentAxes', Axes);
            circle3(cur_click_loc,d,'FaceColor','red');
            last_click_loc=cur_click_loc;
            drawnow
            
            objects_circle_midpoint = [objects_circle_midpoint,cur_click_loc'];
        end
    end
    
    objects_cirlce_radius = repmat(d/2,1,size(objects_circle_midpoint,2));
    
    
    set(Fig,'WindowButtonDownFcn',[]);
    set(Fig,'WindowKeyPressFcn',  []);
    
        function MouseClickCallback(source,eventData)
           cur_click_loc = Axes.CurrentPoint(1,1:2);
        end

        function KeyPressCallback(source,eventData)
            % determine the key that was pressed
            keyPressed = eventData.Key;
            if strcmpi(keyPressed,'escape')
                stop_exec = 1;
            elseif strcmpi(keyPressed,'return')
                not_done = 0;
            end
        end

end

%%
