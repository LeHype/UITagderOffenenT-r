function varargout = set_limits_perc(quant,perc_tb,xyz_axis,LogarithmicPlot)
    arguments
    % mandatory
        quant   double
        perc_tb (1,2) double
    % optional
        xyz_axis        {mustBeMember(xyz_axis,{'x','y','z'})} = 'y'
        LogarithmicPlot {mustBeMember(LogarithmicPlot,[0,1])}  = 0
    end
    switch xyz_axis
        case 'x'
            xyzlim_handle   = @(limits) xlim(limits);
        case 'y'
            xyzlim_handle   = @(limits) ylim(limits);
        case 'z'
            xyzlim_handle   = @(limits) zlim(limits);
    end


    perc_tb = perc_tb/100;
    min_quant = min(reshape(quant,[],1));
    max_quant = max(reshape(quant,[],1));
    if LogarithmicPlot==1
        min_quant = log(min_quant);
        max_quant = log(max_quant);
    end
    delta_quant = max_quant-min_quant;
    lb = min_quant-perc_tb(1)*delta_quant;
    ub = max_quant+perc_tb(2)*delta_quant;
    if LogarithmicPlot==1
        lb = exp(lb);
        ub = exp(ub);
    end
    xyzlim_handle([lb,ub])

    lb_ub = [lb,ub];
    
    if nargout==1
        varargout{1} = lb_ub;
    end

end