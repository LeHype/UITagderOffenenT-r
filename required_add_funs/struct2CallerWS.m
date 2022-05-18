function struct2CallerWS(input_struct,struct_fieldnames)      
    if isa(input_struct,'struct')
        if nargin<2
            struct_fieldnames = fieldnames(input_struct);
        end
        for ii=1:length(struct_fieldnames)
            assignin('caller',struct_fieldnames{ii},input_struct.(struct_fieldnames{ii}));
        end
    else
        warning('Not a struct.')
    end
end