function [rhs_function_handle ,n,m] = getCarDynamics
% model parameter 
l = 1.1;

% model equations
% state_names = {'x','y','v','psi','delta'};
% input_names = {'a','omega'};
rhs_function_handle = @(t,x,u) [ x(3)*cos(x(4))
                                 x(3)*sin(x(4))
                                 u(1)
                                 ...x(3)/l*sin(x(5))
                                 x(3)/l*tan(x(5))
                                 u(2)             ];
n=5;
m=2;
      

end