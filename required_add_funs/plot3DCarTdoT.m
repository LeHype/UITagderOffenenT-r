function car_patch = plotnew3DCar(coords,dims,x,phi,bodyColor,tireColor,plotOuterApprCylinder)

arguments
    %mandatory
        coords
        dims
        x
        phi
    % optional
        bodyColor = 1/250*[0,72,119]  % uds_blue
        tireColor = 1/250*[200,34,84] % uds_red
        plotOuterApprCylinder = 1
end

pos_transl = [x(1:2); 0];
psi_correction = -90*pi/180;
R_psi = rot3D(psi_correction+x(4),'z');
R_delta = rot3D(x(5),'z');
R_phi   = rot3D(phi,'x');


coords.tire_fl = R_delta*R_phi*coords.tire + dims.tire.shift_fl;
coords.tire_fr = R_delta*R_phi*coords.tire + dims.tire.shift_fr;
coords.tire_bl =         R_phi*coords.tire + dims.tire.shift_bl;
coords.tire_br =         R_phi*coords.tire + dims.tire.shift_br;

car_patch(1) = PatchPrism(coords.front,'FaceColor',bodyColor,'Rotation',R_psi,'Translation',pos_transl);
car_patch(2) = PatchCuboid(coords.body,'FaceColor',bodyColor,'Rotation',R_psi,'Translation',pos_transl);

car_patch(3) = PatchCylinder(coords.tire_fl,'Translation',pos_transl,'Rotation',R_psi,'FaceColor',tireColor);
car_patch(4) = PatchCylinder(coords.tire_fr,'Translation',pos_transl,'Rotation',R_psi,'FaceColor',tireColor);
car_patch(5) = PatchCylinder(coords.tire_bl,'Translation',pos_transl,'Rotation',R_psi,'FaceColor',tireColor);
car_patch(6) = PatchCylinder(coords.tire_br,'Translation',pos_transl,'Rotation',R_psi,'FaceColor',tireColor);

car_patch(7) = PatchCylinder(coords.axis_f,'Translation',pos_transl,'Rotation',R_psi,'FaceColor',tireColor);
car_patch(8) = PatchCylinder(coords.axis_b,'Translation',pos_transl,'Rotation',R_psi,'FaceColor',tireColor);



if plotOuterApprCylinder==1
    middleDot = [x(1);x(2);0];
    car_patch(9) = PatchCylinder(coords.YellowCylinder,'Translation',middleDot,'RotCenterShift',[0;0;-0.25],'FaceColor','yellow');
    car_patch(10) = PatchCylinder(coords.OutApprCylinder{1},'Rotation',R_psi,'Translation',pos_transl,'FaceColor',0.4*[1 1 1],'EdgeColor','none');
    car_patch(11) = PatchCylinder(coords.OutApprCylinder{2},'Rotation',R_psi,'Translation',pos_transl,'FaceColor',0.4*[1 1 1],'EdgeColor','none');   
end

end