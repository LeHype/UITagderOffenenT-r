function varnameString = var2varnameString(varargin)
    varnameString = '';
    for ii=1:length(varargin)
        varname = inputname(ii);
        varnameString  = [varnameString,varname,' '];
    end
end