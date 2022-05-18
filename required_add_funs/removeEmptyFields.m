function output_struct = removeEmptyFields(input_struct)
    fn = fieldnames(input_struct);
    for ii=1:length(fn)
        if isempty(input_struct.(fn{ii}))
            input_struct = rmfield(input_struct,fn{ii});
        end
    end
    output_struct = input_struct;
end