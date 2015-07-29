function [x,y] = cnv_thph2xy(th,ph,p)

% Definition of camera parameter vector.
    dx = p(1);
    dy = p(2);
    psi0 = p(3);
    th0 = p(4);
    ph0 = p(5);
    a1 = p(6);
    a2 = p(7)/250; % for scaling
    a3 = p(8)/250^2; % for scaling
    flip = p(9);
    
% Rotate imager-relative zenith and "azimuth" to real zenith and "azimuth"
% angles. ("azimuth" is in quotes because here it is defined as
% north-of-east)
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

% Flip according to back/front-illuminated flag
    if flip
        ph_i = -ph_i;
    end
     
% Mapping from angle north-of-east to polar angle on CCD
    psi = ph_i + (pi - psi0);
    
% Mapping from imager-relative zenith angle to radial distance on CCD
    rho = zeros(size(th_i));
    for i = 1:length(th_i)
        rhos = roots([a3,a2,a1,-th_i(i)]);
        possible = imag(rhos)==0 & rhos > 0;
        sorted = sort(rhos(possible));
        % Give an error if there is no intersection.
        % TODO: make this a more graceful error so that the solver won't
        % fail in these situations, but rather steer away from them.
        if isempty(sorted)
            error('Improper zenith angle falloff: No solution');
        else
            rho(i) = sorted(1); % The smallest intersection
        end
    end
        
% Convert from polar CCD coordinates to cartesian CCD coordinates
    x = rho.*cos(psi) + dx;
    y = rho.*sin(psi) + dy;

end