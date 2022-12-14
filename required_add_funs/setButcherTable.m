function varargout = setButcherTable(integr_method,des_int_order,integr_method_params)
    arguments
        % mandatory
            integr_method
            des_int_order
        % optional 
            integr_method_params struct = [];
    end

    if ~isa(des_int_order,'struct')
        des_int_order_cache = des_int_order;
        clearvars des_int_order
        des_int_order.(integr_method) = des_int_order_cache;
    end

    if ~isempty(integr_method_params)
        struct2CallerWS(integr_method_params)
    else
        eta_ERK2 = 0.5;
        delta_MP = 0;
    end

    %%
    switch integr_method
        case 'ERK'
            if 1
            switch des_int_order.ERK
                case {'1,1','ExplicitEuler'}
                    butcher.name = 'S1O1_ExplEuler_ERK';
                    A = 0;                            
                    b = 1;
                    c = 0;                    
                case {'2,1','Heun_2O'}
                    butcher.name = 'S2O2_Heun_ERK';
                    A = [0, 0;   
                         1, 0];
                    b = [1/2, 1/2];
                    c = [0, 1];                    
                case {'2,2','S2O2_Generic_ERK'}
%                     eta = eta_ERK2;
%                     A = [0, 0;   
%                          eta, 0];
%                     b = [1-1/(2*eta), 1/(2*eta)];
%                     c = [0, eta];
                    A = [0, 0;   
                         1/2+delta_MP, 0];
                    b = [delta_MP/(1+delta_MP), 1/(1+delta_MP)];
                    c = [0, 1/2+delta_MP];
                    butcher.name = 'S2O2_Generic_ERK';
                case {'2,3','S3_DeltaEMP'}
                    A = [0,   0,   0;
                         1,   0,   0;
                         1/4, 1/4, 0];
                    b = [delta_MP, delta_MP, 1-2*delta_MP];
                    c = [0, 1, 1/2];
                    butcher.name = 'S3O2_DeltaEMP_ERK';
                case {'3,1','Heun_3O'} %Heuns Third Order method
                    A = [0,   0,   0;
                         1/3, 0,   0;
                         0,   2/3, 0];
                    b = [1/4, 0, 3/4];
                    c = [0, 1/3, 2/3];
                    butcher.name = 'S3O3_Heun_ERK';
                case '3,2'
                    A = [0,   0,   0;   
                         1/2, 0,   0;
                         -1   2    0];
                    b = [1/6, 4/6, 1/6];
                    c = [0, 1/2, 1];
                case {'4,1','RK4'} %RK4
                    A = [0,   0,   0, 0; 
                         1/2, 0,   0, 0; 
                         0,   1/2, 0, 0; 
                         0,   0,   1, 0];
                    b = [1/6, 1/3, 1/3, 1/6];
                    c = [0, 1/2, 1/2, 1];
                    butcher.name = 'S4O4_RK4_ERK';
                case {'4,2','Ralston_4O'}
                    A = [0    		 0    		0    		0
                         0.4  		 0    		0  			0
                         0.29697761  0.15875964 0    		0
                         0.21810040  -3.05096516 3.83286476 0];
                    b = [0.17476028   -0.55148066  1.20553560  0.17118478];
                    c = [0 0.4 0.45573725 1];
                    butcher.name = 'S4O4_Ralston_ERK';
                case '5,1'
                    A = [ 0 			zeros(6,1);
                          1/5			0				zeros(5,1);
                          3/40			9/40			0				zeros(4,1);
                          44/45			-56/15			32/9			0				zeros(3,1);
                          9372/6561		-25360/2187		64448/6561		-212/729		0				zeros(2,1);
                          9017/3168		-355/33			46732/5247		49/176			-5103/18656		0				Inf;
                          35/384		0				500/1113		125/192			-2187/6784		11/84			0 	 ];
                    b = [37/378	0	250/621	125/594	0	512/1771];
                    c = [0 1/5 3/10 4/5 8/9 1 1];	
                    butcher.name = 'S7O5_Cash_Karp';
                case '6,1'
                    A = [0,       zeros(1,6)
                         1/3,     0,      zeros(1,5)
                         0,       2/3,    0,      zeros(1,4)
                         1/12,    1/3,    -1/12,  0,       zeros(1,3)
                         25/48,   -55/24, 35/48,  15/8,    0,      zeros(1,2)
                         3/20,    -11/20, -1/8,   1/2,     1/10,   0,     zeros(1,1)
                         -261/20, 33/13,  43/156, -118/39, 32/195, 80/39, 0            ];
                    b =  [13/200, 0, 11/40, 11/40, 4/25, 4/25, 13/200];
                    c =  [0, 1/3, 2/3, 1/3, 5/6, 1/6, 1];
                    butcher.name = 'S7O6';
            end
            butcher.A = A;
            butcher.b = b;
            butcher.c = c;
            end
        case 'IRK'
            if 1
            A_expl = [];
            switch des_int_order.IRK
                case {'1,1','ImplicitEuler'}
                    butcher.name = 'S1O1_ImplEuler_IRK';
                    A = 1;                            
                    b = 1;
                    c = 1;                    
                case '2,1'
                    butcher.name = 'S1O2_Implicit_Midpoint_IRK';
                    A = 1/2;
                    b = 1;
                    c = 1/2;                    
                case '3,1'
                    butcher.name = 'S2O3_IRK';
                    A = [1/2+sqrt(3)/6, 0; -sqrt(3)/3, 1/2+sqrt(3)/6];
                    b = [1/2, 1/2];
                    c = [1/2+sqrt(3)/6, 1/2-sqrt(3)/6];
                case '3,2'
                    butcher.name = 'S2O3_Radau_IA_IRK';
                    A = [1/4, -1/4; 1/4, 5/12];
                    b = [1/4, 3/4];
                    c = [1/2+sqrt(3)/6, 1/2-sqrt(3)/6];
                case '4,1' 
                    butcher.name = 'S4O4_IRK';
                    b = [1/12, 5/12, 5/12, 1/12];
                    A = [0,0,0,0;
                         (11-sqrt(5))/120, (25+sqrt(5))/120   , (25-11*sqrt(5))/120, (-1-sqrt(5))/120;
                         (11+sqrt(5))/120, (25+11*sqrt(5))/120, (25-sqrt(5))/120   , (-1+sqrt(5))/120;
                         b];
                    c = [0, (5-sqrt(5))/10, (5+sqrt(5))/10, 1];
                case '4,2'
                    butcher.name = 'S3O4_Lobatto_IIICs_IRK';
                    b = [1/6, 2/3, 1/6];
                    A = [1/6, -1/3, 1/6
                         1/6, 5/12, -1/12;
                         b];
                    c = [0, 1/2, 1];                    
                case '4,3'
                    A = [zeros(1,3);
                         5/24, 1/3, -1/24
                         1/6,  2/3,  1/6 ];
                    b = A(3,:);
                    c = [0, 1/2, 1];
                otherwise
                    error('Choose order for existing method')
            end
            butcher.A = A;
            butcher.b = b;
            butcher.c = c;
            butcher.A_expl = A_expl;
            end
        case 'SPRK'
            if 1
            A_expl = [];
            switch des_int_order.SPRK
                case '2,1'
                    butcher.name = 'Si1O2_IMP_SPRK';
                    A_q = 1/2;
                    A_p = 1/2;
                    b   = 1;
                    c   = 1/2;
                case '2,2'       
                    butcher.name = 'Si2O2_SPRK';
                    gamma = (2-sqrt(2))/2;
                    b_2   = (1-2*gamma)/(4*gamma);
                    A_q = [0,           0,     0;
                           gamma,       gamma, 0;
                           1-b_2-gamma, b_2,     gamma];
                    A_p = A_q;
                    b = [1-b_2-gamma, b_2,     gamma];
                    c = [0,           2*gamma, 1];  
                case '3,1'         
                    butcher.name = 'Si2O3_SPRK';
                    A_q = [ 1/2+sqrt(3)/6, 0; 
                            -sqrt(3)/3,    1/2+sqrt(3)/6 ];
                    A_p = A_q;
                    b   = [1/2, 1/2];
                    c   = [1/2+sqrt(3)/6, 1/2-sqrt(3)/6];
                case '4,1'
                    butcher.name = 'Si3O4_Gaussian_SPRK';
                    A_q = [5/36-sqrt(15)/90    , 2/9-(2*sqrt(15))/45, 5/36-(2*sqrt(15))/45;
                           5/36+sqrt(15)/24    , 2/9                , 5/36-sqrt(15)/24;
                           5/36+(2*sqrt(15))/45, 2/9+(2*sqrt(15))/45, 5/36+sqrt(15)/90];                  
        
                    A_p = [5/36+sqrt(15)/90    , 2/9-sqrt(15)/15    , 5/36-(2*sqrt(15))/45;
                           5/36+sqrt(15)/36    , 2/9                , 5/36-sqrt(15)/36;
                           5/36+(2*sqrt(15))/45, 2/9+sqrt(15)/15    , 5/36-sqrt(15)/90];
                    b   = [5/18, 4/9, 5/18];
                    c   = [(5-sqrt(15))/10, 1/2, (5+sqrt(15))/10]; 
                case '4,2'
                    butcher.name = 'Si3O4_Lobatto_SPRK';
                    b = [1/12, 5/12, 5/12, 1/12];
                    c = [0, (5-sqrt(5))/10, (5+sqrt(5))/10, 1];
                    A_p = [0,0,0,0;
                           (11-sqrt(5))/120, (25+sqrt(5))/120   , (25-11*sqrt(5))/120, (-1-sqrt(5))/120;
                           (11+sqrt(5))/120, (25+11*sqrt(5))/120, (25-sqrt(5))/120   , (-1+sqrt(5))/120;
                           b];
                    A_q = [1/12, (-1+sqrt(5))/24 ,    (-1-sqrt(5))/24    , 0;
                           1/12, (25-sqrt(5))/120,    (25-11*sqrt(5))/120, 0;
                           1/12, (25+11*sqrt(5))/120, (25+sqrt(5))/120   , 0;
                           1/12, (11+sqrt(5))/24,     (11-sqrt(5))/24    , 0];
                case '6,1'
                    butcher.name = 'Si4O6_SPRK';
                    A_q = [ 1/12, -sqrt(5)/12,       sqrt(5)/12,        -1/12;
                            1/12, 1/4,               (10-7*sqrt(5))/60, sqrt(5)/60;
                            1/12, (10+7*sqrt(5))/60, 1/4,               -sqrt(5)/60;
                            1/12, 5/12,              5/12,              1/12         ];
                    A_p = [ zeros(1,4);
                            (5+sqrt(5))/60, 1/6,               (15-7*sqrt(5))/60, 0;
                            (5-sqrt(5))/60, (15+7*sqrt(5))/60, 1/6,               0;
                            1/6,            (5-sqrt(5))/12,    (5+sqrt(5))/12,    0  ];
                    b = [1/12, 5/12, 5/12, 1/12];
                    c = [0, (5-sqrt(5))/10, (5+sqrt(5))/10, 1];
                case '8,1'
                    butcher.name = 'Si6O8_SPRK';
                    A_q = [ 1/20, -7/60,                   2/15,                 -7/60,                   1/20;
                            1/20, 29/180,                  (47-15*sqrt(21))/315, (203-30*sqrt(21))/1260,  -3/140;
                            1/20, (329+105*sqrt(21))/2880, 73/360,               (329-105*sqrt(21))/2880, 3/160;
                            1/20, (203+30*sqrt(21))/1260,  (47+15*sqrt(21))/315, 29/180,                  -3/140;
                            1/20, 49/180,                  16/45,                49/180,                  1/20    ];
                    b   = [1/20, 49/180, 16/45, 49/180, 1/20];
                    A_p = construct_barA(A_q,b);
                    c   = [0 1/2-sqrt(21)/14, 1/2, 1/2+sqrt(21)/14, 1];
                otherwise
                    error('Choose order for existing method')
            end
            SPRK_check(A_p,b,A_q,b)
            butcher.A_q = A_q;
            butcher.A_p = A_p;
            butcher.b   = b;
            butcher.c   = c;
            butcher.A_expl = A_expl;
            end
        case 'ERK_adaptive'
            if 1
            switch des_int_order.ERK_adaptive
                case '12,1'
                    A = [0, 0;   
                         1, 0];
                    b = [1/2, 1/2];
                    c = [0, 1];
                    b_e = [1 0]-b;
                case '45,1'
                    A_5 = [ zeros(1,7);
                            1/5          0          0			0          0           0     0;
                            3/40         9/40		0           0          0           0     0;
                            44/45       -56/15      32/9		0          0           0     0;
                            19372/6561	-25360/2187	64448/6561	-212/729   0           0     0;
                            9017/3168	-355/33     46732/5247	49/176    -5103/18656  0     0;	
                            35/384       0          500/1113	125/192	  -2187/6784   11/84 0  ];
                    b_5 = [ 5179/57600 0  7571/16695 393/640 -92097/339200 187/2100	  1/40 ];
                    c_5 = [0 1/5 3/10 4/5 8/9 1 1];
            %       A_4 = A_5(1:6,6);
            %       b_4 = [ 35/384     0  500/1113	 125/192 -2187/6784    11/84           ];
            %       c_4 = c_5(1:6);
                    b_e = [ 71/57600   0 -71/16695   71/1920 -17253/339200 22/525    -1/40 ];
                    A = A_5;
                    b = b_5;
                    c = c_5;
            end
            butcher.A = A;
            butcher.b = b;
            butcher.c = c;
            butcher.b_e = b_e;
            end
        case 'IMP'
            butcher.A = 1/2;
            butcher.b = 1;
            butcher.c = 1/2;
    end

    %%
    butcher.s = length(butcher.c);
    varargout{1} = butcher;
    varargout{2} = butcher.s;
end