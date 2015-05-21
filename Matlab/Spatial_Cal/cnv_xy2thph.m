function [th,ph] = cnv_xy2thph(x,y,dx,dy,psi0,th0,ph0,a,flip)

    rho = sqrt((x-dx).^2 + (y-dy).^2);
    psi = atan2(y-dy,x-dx);
    th_i = a*rho;
    ph_i = psi - (pi - psi0);
    if flip
        ph_i = -ph_i;
    end
    r_i = [sin(th_i).*cos(ph_i);
    	   sin(th_i).*sin(ph_i);
    	   cos(th_i)];
    Ry = [cos(th0)  0 sin(th0);
          0         1 0;
          -sin(th0) 0 cos(th0)];
    Rz = [cos(ph0) -sin(ph0) 0;
          sin(ph0) cos(ph0)  0;
          0        0         1];
    r_l = Rz*Ry*r_i;
    ph = atan2(r_l(2,:),r_l(1,:));
    th = acos(r_l(3,:));
    
end