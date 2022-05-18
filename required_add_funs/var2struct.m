function output_struct = var2struct(varargin)
    for ii=1:length(varargin)
        var_name = inputname(ii);
        output_struct.(var_name) = varargin{ii};
    end
end