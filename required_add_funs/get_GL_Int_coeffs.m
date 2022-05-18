function [alpha_vec,x_vec] = get_GL_Int_coeffs(n_GL)
    switch n_GL
        case 1
            x_vec     = 0;
            alpha_vec = 2;
        case 2
            x_vec     = [-sqrt(1/3), sqrt(1/3)]';
            alpha_vec = [1, 1]';      
        case 3
            x_vec     = [-sqrt(3/5), 0, sqrt(3/5)]';
            alpha_vec = [5/9, 8/9, 5/9]';  
        case 4
            x_vec     = [-sqrt(3/7+2/7*sqrt(6/5)), -sqrt(3/7-2/7*sqrt(6/5)), sqrt(3/7-2/7*sqrt(6/5)), sqrt(3/7+2/7*sqrt(6/5))]';
            alpha_vec = [(18-sqrt(30))/36, (18+sqrt(30))/36, (18+sqrt(30))/36, (18-sqrt(30))/36];  
%             case 5
%                 c_vec     = 0;
%                 alpha_vec = 2;  
        otherwise
            error('Choose n smaller than 5')
    end
end