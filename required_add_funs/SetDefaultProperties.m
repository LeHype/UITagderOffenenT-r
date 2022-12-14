function SetDefaultProperties(options)
    arguments
        % mandatory
            % nothing
        % optional 
            options.Fontsize         = 15;
            options.Linewidth        = 0.5; %default
            options.Interpreter      {mustBeMember(options.Interpreter,{'tex','latex','none'})} = 'latex'
            options.FactorySettings  {mustBeMember(options.FactorySettings,[0,1])} = 0
    end
    struct2CallerWS(options)
    
    groot_objects = repmat({groot},1,3);
    zero_objects  = repmat({0},1,7);
    objects  = [groot_objects,zero_objects];
    NV_pairs = {'AxesTitleFontSizeMultiplier', 1;
                'AxesLabelFontSizeMultiplier', 1;
                'AxesTickLabelInterpreter',    Interpreter; 
                'TextInterpreter',             Interpreter;
                'LegendInterpreter',           Interpreter;
                'AxesFontSize',                Fontsize;
                'LegendFontSize',              Fontsize;
%                 'AxesLineWidth',               Linewidth;
                'LineLineWidth',               Linewidth;
                'LegendFontSizeMode',          'auto';%'manual';
                'LegendAutoUpdate',            'off'};
 
    for ii=1:length(objects)
        if FactorySettings==0
            Value_ii = NV_pairs{ii,2};
        else
            Value_ii = get(objects{ii},['Factory',NV_pairs{ii,1}]);
        end
        set(objects{ii},['Default',NV_pairs{ii,1}],Value_ii)
    end

end