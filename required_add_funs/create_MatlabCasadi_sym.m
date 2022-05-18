function create_MatlabCasadi_sym(MatOrCas,var_texts,workspace_name,Casadi_var_type)
    % example calls:
    %   create_MatlabCasadi_sym('Matlab',{'t','h','x_0 5 1','x_1 5 1','u_0','u_1'},'caller')
    %   create_MatlabCasadi_sym('Casadi',{'t','h','x_0 5 1','x_1 5 1','u_0','u_1'},'MX')
    % if variable is scalar you can skip dimension specification
    var_names = var_texts;
    var_dims  = var_texts;
    for ii=1:length(var_names)
        var_text_ii = var_names{ii};
        ind = strfind(var_text_ii,' ');
        if isempty(ind)
            var_names(ii) = {var_text_ii };
            var_dims(ii)  = {'1 1'};
        else
            ind = ind(1);
            var_names(ii) = {var_text_ii(1:ind-1)};
            var_dims(ii)  = {var_text_ii(ind+1:end)};
        end
    end
    Matlab_var = strcmp(MatOrCas,'Matlab');
    Casadi_var = strcmp(MatOrCas,'Casadi');
    code_real  = ',''real''';
    code_3     = ');';
    for ii=1:length(var_names)
        code_1 = [var_names{ii},' = '];
        code_2 = ['sym(''',var_names{ii},''',[',var_dims{ii},']'];
        if Matlab_var==1
            if strcmp(var_dims{ii},'1 1')
                evalc(['syms ',var_names{ii},' real']);
            else
                evalc([code_1,code_2,code_real,code_3]);
            end
        elseif Casadi_var==1
            if nargin<=3
                Casadi_var_type = 'SX';
            elseif strcmp(Casadi_var_type,{'MX','SX'})~=1
                error([Casadi_var_type,' is not a valid Casadi variable type.'])
            end
            code = [code_1,Casadi_var_type,'.',code_2,code_3];
            try
                evalc(code);
            catch
                import casadi.*
                evalc(code);
            end
        else
            error('Choose between "Matlab" and "Casadi".')
        end
        if nargin<=2 || isempty(workspace_name)
            workspace_name = 'base';
        end
            
        % assign variables to base workspace
            evalc(['assignin(''',workspace_name,''',''',var_names{ii},''',',var_names{ii},');']);
    end
end