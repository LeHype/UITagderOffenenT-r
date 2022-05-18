function lb_ub = set_limits_perc(quant,perc_tb,xyzlim_handle,LogarithmicPlot)
    arguments
    % mandatory
        quant   double
        perc_tb (1,2) double
    % optional
        xyzlim_handle   = @(limits) ylim(limits);
        LogarithmicPlot = 0;
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
end