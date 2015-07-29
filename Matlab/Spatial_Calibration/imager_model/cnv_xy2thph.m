function [th,ph] = cnv_xy2thph(x,y,p)

% Definition of camera parameter vector p
    dx = p(1);
    dy = p(2);
    psi0 = p(3);
    th0 = p(4);
    ph0 = p(5);
    a1 = p(6);
    a2 = p(7)/250; % for scaling
    a3 = p(8)/250^2; % for scaling
    flip = p(9);
    
% Convert cartesian CCD coordinates to polar CCD coordinates
    rho = sqrt((x-dx).^2 + (y-dy).^2);
    psi = atan2(y-dy,x-dx);

% Mapping from radial distance on CCD to imager-relative zenith angle
    th_i = a1*rho + a2*rho.^2 + a3*rho.^3;

% Mapping from polar angle on CCD to angle north-of-east
    ph_i = psi - (pi - psi0);
    
% Flip according to back/front-illuminated flag
    if flip
        ph_i = -ph_i;
    end
    
% Rotate imager-relative zenith and "azimuth" to real zenith and "azimuth"
% angles. ("azimuth" is in quotes because here it is defined as
% north-of-east)
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