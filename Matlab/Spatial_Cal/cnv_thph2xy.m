function [x,y] = cnv_thph2xy(th,ph,dx,dy,psi0,th0,ph0,a,flip)
    
    r_l = [sin(th).*cos(ph);
    	   sin(th).*sin(ph);
    	   cos(th)];
    Ry = [cos(th0)  0 -sin(th0);
          0         1 0;
          sin(th0) 0 cos(th0)];
    Rz = [ cos(ph0) sin(ph0) 0;
          -sin(ph0) cos(ph0)  0;
          0        0         1];  
    r_i = Ry*Rz*r_l;
    th_i = acos(r_i(3,:));
    ph_i = atan2(r_i(2,:),r_i(1,:));
    if flip
        ph_i = -ph_i;
    end
    psi = ph_i + (pi - psi0);
    rho = th_i/a;
    x = rho.*cos(psi) + dx;
    y = rho.*sin(psi) + dy;

end