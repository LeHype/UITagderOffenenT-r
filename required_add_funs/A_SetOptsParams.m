function WS_struct = A_SetOptsParams(OptsParams)
    arguments
        % mandatory
            %nothing
        % optional 
            OptsParams.vars_as_block   {mustBeMember(OptsParams.vars_as_block,[0,1])}   = 0; 
            OptsParams.constr_X_interm {mustBeMember(OptsParams.constr_X_interm,[0,1])} = 0;
            OptsParams.scale_x         {mustBeMember(OptsParams.scale_x,[0,1])}         = 0;
            OptsParams.scale_u         {mustBeMember(OptsParams.scale_u,[0,1])}         = 0;
            OptsParams.parametrize_OP  {mustBeMember(OptsParams.parametrize_OP,[0,1])}  = 1;
            OptsParams.option_set      {mustBeMember(OptsParams.option_set,1:3)}        = 1;
            OptsParams.t_grid_normed_opt {mustBeMember(OptsParams.t_grid_normed_opt,1:3)} = 1;
            OptsParams.casadi_options struct = [];  

            % solver options                       
            OptsParams.solver_choice  (1,:) char = 'ipopt'; %'ipopt','worhp'
            OptsParams.solver_options struct = []; 

            % integrator choice + cost functional
            OptsParams.integr_method {mustBeMember(OptsParams.integr_method,{'IRK','IMP','ERK','ERK_adaptive','SPRK'})} = 'ERK';              
            OptsParams.L_int_type    {mustBeMember(OptsParams.L_int_type,{'Mayer','GL','analytic'})} = 'Mayer';  % analytic does not consider x_ki yet!
            OptsParams.deltaEMP_Mayer_on  = 0;

            % integrator choice                                
            OptsParams.des_int_order   struct = struct('SPRK','3,1','IRK','2,1','ERK','4,1');            
            OptsParams.integr_method_params struct = [];   
            OptsParams.n_GL {mustBeMember(OptsParams.n_GL,[1 2 3 4])} = 2; 
            
            % PATH FOLLOWING
                OptsParams.PathFollowing {mustBeMember(OptsParams.PathFollowing,[0,1])} = 0; 
    end
    struct2CallerWS(OptsParams)

    if isempty(casadi_options)
        casadi_options = struct;
    end
    if ~isempty(solver_options)
        casadi_options.(solver_choice) = solver_options;
    end

    %% problem formulation                              
    switch option_set
        case 1 % classical ZOH
            reform_dyn_xu = 0;
            OC_U_ZOH      = 1;
            simCT_U_FOH   = 0; 
        case 2 % consider FOH explicitely in OCP
            reform_dyn_xu = 0;
            OC_U_ZOH      = 0;
            simCT_U_FOH   = 1; %if OC_U_ZOH=1 and simCT_U_FOH=1 the computed inputs are connected linearly
        case 3 % augment dynamics and consider ZOH in OCP, then connect computed inputs (u, not du) linearly
            reform_dyn_xu = 1;
            OC_U_ZOH      = 1;
            simCT_U_FOH   = 1; %if OC_U_ZOH=1 and simCT_U_FOH=1 the computed inputs are connected linearly
    end   
    

    %% output WS_struct
    clearvars OptsParams ans
    strWS2WS_struct = who();
    strWS2WS_struct = strjoin(strWS2WS_struct,',');
    WS_struct = eval(sprintf('var2struct(%s);',strWS2WS_struct));
    try 
        WS_struct = rmfield(WS_struct,'ans');
    end
end