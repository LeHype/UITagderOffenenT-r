function [dims,coords] = createNew3DCarObject(BodyLength,absb_delta,options)
arguments
    %mandatory
        BodyLength
        absb_delta
    % optional
        options.BodyWidth   = 1
        options.BodyHeight  = 1.05
        options.FrontLength = 0.5
        options.TireWidth   = 0.1
        options.TireRadius  = 0.2
end
struct2CallerWS(options)


%% dimensions
dims.body.l = BodyWidth;
dims.body.d = BodyLength;
dims.body.h = BodyHeight;

dims.front.h = BodyHeight;
dims.front.g = BodyWidth;
dims.front.d = FrontLength;

dims.tire.r = TireRadius;
dims.tire.h = TireWidth;
dims.tire.N = 11;

dims.dh = dims.tire.r - 0.5*dims.tire.r;
dims.tire_gap = -TireWidth/2;%dims.tire.h;

dims.axis.r = 0.1;
dims.axis.h = dims.body.l+2*(dims.tire_gap+dims.tire.h);
dims.axis.N = 11;
    
dims.x_tire = (dims.body.l+dims.tire.h+2*dims.tire_gap)/2;
shiftTireFront = BodyLength*1/5;
shiftTireBack  = BodyLength*4/5;


tg = dims.front.g;
th = dims.front.d;
td = dims.front.h;

ldh_outer = [dims.axis.h; 
             dims.front.d+dims.body.d-2*dims.axis.r+dims.tire.r;
             dims.dh+dims.body.h];

R_delta_OutAppr = rot3D(absb_delta,'z');

% coords
    coords.front = CoordsPrism(tg,td,th,'RotCenterShift',[0 0*td/2 0]','Rotation',rot3D(-90,'x','Unit','deg'), ...
        'Translation',[0 shiftTireFront dims.front.h/2+dims.dh]');
    %coords.front = PatchPrism(dims.front,'RotCenterShift',[0 0*td/2 0]','Rotation',rot3D(-90,'x','Unit','deg'), ...
     %   'Translation',[0 0 dims.front.h/2+dims.dh]');
    coords.body = CoordsCuboid(dims.body.l,dims.body.d,dims.body.h, ...
        'Translation',[0 shiftTireFront-dims.body.d/2 dims.body.h/2+dims.dh]');
    coords.tire = CoordsCylinder(dims.tire.r,dims.tire.h,dims.tire.N,'Rotation',rot3D(-90,'y','Unit','deg'), ...
        'Translation',[0 0*shiftTireFront 0]');%,[0 0 dims.tire.r]');
    coords.axis = CoordsCylinder(dims.axis.r,dims.axis.h,dims.axis.N,'Rotation',rot3D(-90,'y','Unit','deg'), ...
        'Translation',[0 shiftTireFront dims.tire.r]');
    coords.outer = CoordsCuboid(ldh_outer(1),ldh_outer(2),ldh_outer(3), ...
        'Translation',[0 shiftTireFront-ldh_outer(2)/2+dims.front.d ldh_outer(3)/2]');

    coords.OutApprCylinder(1) = {CoordsCylinder(1.5/2,BodyHeight,31,'Translation',[0; -FrontLength - BodyLength*1.2/10 + (FrontLength+BodyLength*1/5); BodyHeight/2])};
    coords.OutApprCylinder(2) = {CoordsCylinder(1.5/2,BodyHeight,31,'Translation',[0; -FrontLength - BodyLength*7.2/10 + (FrontLength+BodyLength*1/5); BodyHeight/2])};
    % circle_top    = [0;-FrontLength - BodyLength*1.2/10];
    % circle_bottom = [0;-FrontLength - BodyLength*7.2/10];
    
    dims.tire.shift_fl = [-dims.x_tire -0*shiftTireFront dims.tire.r]';
    dims.tire.shift_fr = [+dims.x_tire -0*shiftTireFront dims.tire.r]';
    dims.tire.shift_bl = [-dims.x_tire -shiftTireBack+shiftTireFront dims.tire.r]';
    dims.tire.shift_br = [+dims.x_tire -shiftTireBack+shiftTireFront dims.tire.r]';

    coords.tire_fl = R_delta_OutAppr*coords.tire  + dims.tire.shift_fl;
    coords.tire_fr = R_delta_OutAppr'*coords.tire + dims.tire.shift_fr;
    coords.tire_bl =                  coords.tire + dims.tire.shift_bl;
    coords.tire_br =                  coords.tire + dims.tire.shift_br;
    
    coords.axis_f = coords.axis + [0 -shiftTireFront 0]';
    coords.axis_b = coords.axis + [0 -shiftTireBack  0]';

    

end