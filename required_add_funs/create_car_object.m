function [car,car_colors,indiv_part,circles,body_center] = create_car_object(main_color)
%% create car object
l_front = 0.5;
w_body = 1;
l_body = 1.8;
w_tire = 0.1;
l_tire = 0.4;
circle_diam = [1.5, 1.5];

front = [ -w_body/2,  w_body/2, 0
          -l_front,  -l_front,  0 ];
body = [ -w_body/2,        w_body/2,        w_body/2, -w_body/2
         -l_front-l_body, -l_front-l_body, -l_front,  -l_front  ];
tire_0 = [ -w_tire/2,  w_tire/2, w_tire/2, -w_tire/2 
           -l_tire/2, -l_tire/2, l_tire/2,  l_tire/2 ];
tire_fr = [ (w_body+0*w_tire)/2; -(l_front+l_body*1/5)] + tire_0;
tire_br = [ (w_body+0*w_tire)/2; -(l_front+l_body*4/5)] + tire_0;
tire_fl = [-(w_body+0*w_tire)/2; -(l_front+l_body*1/5)] + tire_0;
tire_bl = [-(w_body+0*w_tire)/2; -(l_front+l_body*4/5)] + tire_0;

circle_top    = [0;-l_front - l_body*1.2/10];
circle_bottom = [0;-l_front - l_body*7.2/10];

car = {front,body,tire_fr,tire_br,tire_fl,tire_bl};
uds_blue = 1/250*[0,72,119];
uds_red  = 1/250*[200,34,84];
if nargin==0
    main_color = 'blue';
end
switch main_color 
    case 'blue'
        main_color_apply = uds_blue;
    case 'red'
        main_color_apply = uds_red;
    otherwise
        main_color_apply = uds_blue;
end
% car_colors = {'blue','blue','black','black','black','black'};
car_colors = {main_color_apply,main_color_apply,'black','black','black','black'};


%% move midpoint between front tires to (0,0)
car           = rot_transl_object(car,0,'deg',        [0;(l_front+l_body*1/5)],[],[]);
circle_top    = rot_transl_part(circle_top,0,'deg',   [0;(l_front+l_body*1/5)]);
circle_bottom = rot_transl_part(circle_bottom,0,'deg',[0;(l_front+l_body*1/5)]);

tire_fr = car{3};
tire_fl = car{5};
tire_fr_M = [(tire_fr(1,1)+tire_fr(1,2))/2; (tire_fr(2,2)+tire_fr(2,3))/2];
tire_fl_M = [(tire_fl(1,1)+tire_fl(1,2))/2; (tire_fl(2,2)+tire_fl(2,3))/2];

indiv_part.nr    = [3,5];
indiv_part.angle = [0,0];
indiv_part.axis  = [tire_fr_M,tire_fl_M];

circles.midpoint = [circle_top, circle_bottom];
circles.radius   = 1/2*circle_diam;

body_center = [0 0]';

end