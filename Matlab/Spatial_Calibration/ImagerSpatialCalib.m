function [ccd_az, ccd_el, cost, p] = ...
        ImagerSpatialCalib(d, lat, lon, year, doy, time, p0, ... 
                 numAvailableStars, verbose, useold, varargin)
% [ccd_az, ccd_el, cost, p] = ...
%        ImagerSpatialCalib(d, lat, lon, year, doy, time, p0, ... 
%                 numAvailableStars, verbose, useold [, known])
%
% Match the stars in the image d to a known star field in order to infer
% the characteristics of the imager.  Use these characteristics to return
% the elevation and azimuth mappings.
% v1.00: Original version
% v1.01: Go from linear falloff to cubic polynomial.
%        Include "useold" parameter.
% v1.02: Allow user to specify previously known camera parameters.
%        Always return ccd_az in the range [0,360).
%
% Inputs:
% d - MxN sky image
% lat - latitude of imager
% lon - longitude of imager
% year - the year the image was taken (UT)
% doy - the day of year the image was taken (1-366) (UT)
% time - the time of day (in seconds) the image was taken (UT)
% p0 - (9x1) The initial guess for the camera parameter vector (see below
%       for camera parameter vector definition).
% numAvailableStars - the number of stars in the catalog to use. (Typical:
%       100 for allsky, 1000 for narrow-field). Instead of a number, the
%       string 'all' is also accepted, in which case all stars in the 
%       catalog will be shown.
% verbose - If true, show the image, stars, and star field for visual
%           confirmation after convergence.
% useold - If true, don't prompt the user to pick stars, but rather use the
%          last-picked stars.  This is useful for quick debugging.
% known - (optional) (9x1) A logical vector that indicates which elements
%         of the initial guess p0 are known a priori, meaning they don't
%         have to be optimized.  It is always assumed that known(end) =
%         True, since the flip flag has to be correct in order to visually
%         match up stars.
%
% Outputs:
% ccd_az - MxM mapping from pixel to azimuth, in degrees
% ccd_el - MxM mapping from pixel to elevation, in degrees
% cost - A scalar measure of how "far" the image stars are from the
%        known star field, in pixels.  Lower is better.
% p - The solved imager parameter vector: [dx,dy,psi0,th0,ph0,a1,a2,a3,flip]
%     (See below). This is probably not useful to the user.
%
% Definition of camera parameter vector:
%           p = [dx,dy,psi0,th0,ph0,a1,a2,a3,flip]
% (dx,dy) - pix - the pixel that lies on the optical axis (fractional 
%    values allowed)
% psi0 - rad - the azimuthal direction on the CCD that corresponds to the
%   "minus zenith angle" direction in real space.  The angle definition
%   is psi0 = atan2(y-dy,x-dx) where [x-dx,y-dy] is a vector that points
%   "up" from the optical axis on the CCD. This is a rotation of the
%   imager about the optical axis.  Note that when using the function
%   "image", the Y direction is displayed in reverse, but the math is
%   unaffected. Note that for all-sky imagers, there is a degree of
%   freedom between psi0 and ph0.
% th0 - rad - the zenith angle look direction of the optical axis
% ph0 - rad - the azimuthal angle look direction of the optical axis. Note
%   that the coordinate convention used here is that East is +x and North
%   is +y. This makes th0 and ph0 consistent with the typical math
%   definitions of spherical coordinates.
% a1,a2,a3 - the coefficients of the cubic mapping from radial distance 
%   on the CCD to zenith angle.
%   These parameters describe the mapping from real space to CCD
%   coordinates.
%   zenith = a1*r + a2*r^2/250 + a3*r^3/250^2
% flip - a binary (0 or 1) flag describing the handedness of the
%   imager, which corresponds in some way to whether the CCD is back- or
%   front-illuminated.
%
% Algorithm: Display 2 pictures: one is the image, and one is where the
% stars are expected to be from the initial guess p0. The size of circles
% in the second image correspond to the expected brightnesses of the stars.
% The user clicks on some stars in the image, presses enter, and then
% clicks on the corresponding stars in the star field.  From this, the
% optimal camera parameters are calculated and elevation/azimuth mappings
% returned.
%
% Usage:
%   The first 6 parameters in the initial guess do not need to be close to
% the true values for the algorithm to converge, but they probably need to
% be close enough in order to manually select the correct stars.  The 7th
% parameter (the flip flag) needs to be correct for the algorithm to
% converge.
%   If stars are hard to differentiate from noise, try opening two
% images and flipping back and forth.
%
% Author: Brian Harding, University of Illinois at Urbana-Champaign, 2012

% Future Work:
% 1) Automate image star selection.  This can be severely affected by
% noise.
% 2) Automate the inversion by removing manual catalog star selection. This
% will require constructing an objective function that is convex (which I
% think is impossible), or coming up with some heuristic to search
% the parameter space for the global minimum (which I've found is
% unreliable).

% Roadmap:
% (1) From star catalog, find elev/azim of the brightest stars.
% (2) Open a figure that has the image, and the catalog stars.
% (3) Let the user click on the stars to match up.
% (4) Find the set of camera parameters p that optimizes star match.
% (5) Calculate and return the azimuth and elevation mappings from p.
% (6) If requested, plot the star fields

% Settings
    plotprc = [40 99]; % percentiles that define the max and min of image

% Parse inputs
    known = false(size(p0));
    if (~isempty(varargin))
        known = varargin{1};
    end
    known(end) = true; % we have to assume the flip flag is correct
        

% (1) From star catalog, find elev/azim of the brightest stars.
    N = numAvailableStars; % Use the N brightest stars in the catalog
    el_cutoff = 0; % Stars in the catalog below this elevation angle will 
                    % not be used
    cat = load('yale_catalogue');
    [cat_az_all,cat_el_all] = cnv_cel2azel(year,doy,time,cat.obj_ra,...
        cat.obj_dec, lat, lon);
    m_el = find(cat_el_all > el_cutoff);
    [~,m_br] = sort(cat.obj_mag(m_el));
    m = m_el(m_br); % "m" indexes cat_az and cat_el
    cat_az = cat_az_all(m);
    cat_el = cat_el_all(m);
    cat_br = cat.obj_mag(m);
    % Take N brightest
    if ~strcmp(N,'all')
        cat_az = cat_az(1:N);
        cat_el = cat_el(1:N);
        cat_br = cat_br(1:N);
    end
    % Convert to radians
    cat_az = (pi/180)*cat_az;
    cat_el = (pi/180)*cat_el;
    cat_az_all = (pi/180)*cat_az_all;
    cat_el_all = (pi/180)*cat_el_all;

    
% (2) Open a figure that has the image and the catalog stars
    if (~useold)
        h = figure;
        subplot(1,2,1);
        imagesc(d,prctile(d(:),plotprc));
        colormap('gray'); 
        axis equal
        axis tight;
        hold on;
        for r = [pi/6, pi/3, pi/2];
            tmp_ph = linspace(0,2*pi);
            tmp_th = 0*tmp_ph + r;
            [xxx,yyy] = cnv_thph2xy(tmp_th,tmp_ph,p0);
            plot(xxx,yyy,'g-');
        end


        subplot(1,2,2);
        cat_th = pi/2 - cat_el;
        cat_ph = pi/2 - cat_az;
        [xx,yy] = cnv_thph2xy(cat_th,cat_ph,p0);
        s = 2.512.^(5.5 - cat_br);
        scatter(xx,yy,s,'k');
        set(gca,'YDir','Reverse');
        hold on;
        % plot circles
        for r = [pi/6, pi/3, pi/2];
            tmp_ph = linspace(0,2*pi);
            tmp_th = 0*tmp_ph + r;
            [xxx,yyy] = cnv_thph2xy(tmp_th,tmp_ph,p0);
            plot(xxx,yyy,'k-');
        end
        % Note North and East
        th = 70*pi/180;
        ph = pi/2; % north
        [xxx,yyy] = cnv_thph2xy(th,ph,p0);
        text(xxx,yyy,'N','fontsize',15,'color','k','HorizontalAlignment',...
            'center','VerticalAlignment','middle')
        ph = 0; % north
        [xxx,yyy] = cnv_thph2xy(th,ph,p0);
        text(xxx,yyy,'E','fontsize',15,'color','k','HorizontalAlignment',...
            'center','VerticalAlignment','middle')
        
        hold off;
        axis equal
        axis tight
        xlim([0 size(d,2)]);
        ylim([0 size(d,1)]);
        view(0,90);
    
% (3) Let the user click on the stars to match up.
        subplot(1,2,1);
        title(sprintf('Click on stars, then press enter.'));
        starxy = [];
        kk = 1;
        xy = ginput(1);
        hold on;
        while (~isempty(xy))
            x = xy(1);
            y = xy(2);
            %         % look around a little, so the user doesn't have to 
            %         % click right on the star
            %         xy = round(xy);
            %         rr = 2; % radius of box to search around
            %         del = -rr:rr;
            %         window = d(xy(2)+del, xy(1)+del);
            %         [~,idx] = max(window(:));
            %         [dy,dx] = ind2sub(size(window),idx);
            %         y = xy(2) + dy - rr - 1;
            %         x = xy(1) + dx - rr - 1;
            % record this star and plot a marker
            starxy(kk,:) = [x y];
            plot(x,y,'go');
            text(x,y,num2str(kk),'color',[1 1 1]);
            kk = kk + 1;
            xy = ginput(1);
        end
        hold off;
        star_x = starxy(:,1)';
        star_y = starxy(:,2)';
        L = length(star_x);
        title ''

        subplot(1,2,2);
        title(sprintf('Click on %i corresponding stars (in order)',L));
        cat_idx = [];
        cat_orig_idx = [];
        hold on;
        for kk = 1:L
            xy = ginput(1);
            x = xy(1);
            y = xy(2);
            % find cat star closest to the selected point
            [~,i] = min((xx-x).^2 + (yy-y).^2);
            scatter(xx(i),yy(i),s(i),'go');
            cat_idx(kk) = i;
            text(x,y,num2str(kk));
            % also calculate the original index for later retrieval
            [~,origi] = min(abs(cat_az_all - cat_az(i)) + ...
                            abs(cat_el_all - cat_el(i)));
            cat_orig_idx(kk) = origi;
        end
        cat_azi = cat_az(cat_idx);
        cat_eli = cat_el(cat_idx);

        close(h);
        pause(0.3); % to let the figure close

        % Save the picked stars to possible use in the next call
        assignin('base','cat_orig_idx',cat_orig_idx);
        assignin('base','star_x',star_x);
        assignin('base','star_y',star_y);

    else % useold is True
        try
            cat_orig_idx = evalin('base','cat_orig_idx');
            star_x = evalin('base','star_x');
            star_y = evalin('base','star_y');
        catch e
            error('No saved stars: set useold==0');
        end
            
        cat_azi = cat_az_all(cat_orig_idx);
        cat_eli = cat_el_all(cat_orig_idx);

    end

% (4) Find the set of camera parameters p that optimizes star match.
    % Use parameter scaling, and don't solve for known parameters
    scale = [100 100 1 1 1 0.01 0.001 0.001 1]'; % scaling for numerical issues
    pguess = p0./scale;
    f = @(p_est) star_fit_resid(scale.*replace(pguess,p_est,~known),...
                            cat_azi,cat_eli,star_x,star_y)';
%     [p1,cost] = LMFnlsq(f,pguess(~known),'MaxIter',50000);
    opts = optimset('lsqnonlin');
    opts.MaxFunEvals = 50000;
    opts.MaxIter = 5000;
    [p1,cost] = lsqnonlin(f,pguess(~known),[],[],opts);
    p = pguess;
    p(~known) = p1;
    p = p.*scale;
    
% (5) Calculate and return the azimuth and elevation mappings from p.
    [X,Y] = meshgrid(1:size(d,2),1:size(d,1));
    [th,ph] = cnv_xy2thph(X(:)',Y(:)',p);
    ccd_el = zeros(size(d));
    ccd_az = zeros(size(d));
    ccd_el(:) = 90 - 180/pi * th;
    ccd_az(:) = 90 - 180/pi * ph;
    ccd_az = mod(ccd_az,360);

% (6) If requested, plot the star fields
    if verbose
        figure;
        imagesc(d,prctile(d(:),plotprc)); %scale image
        colormap gray
        hold on;
        plot(star_x,star_y,'go','MarkerSize',8);
        % TODO: plot more cat stars
        [cat_x,cat_y] = cnv_thph2xy(pi/2-cat_eli, pi/2-cat_azi,p);
        plot(cat_x,cat_y,'m^','MarkerSize',5);
        hold off;
    end

end



function f = star_fit_resid(p,cat_az,cat_el,star_x,star_y)
% Return a vector of distances that represent the "badness of fit" between
% the stars in the catalog (given in azimuth/elevation in radians) and the
% stars in the image (given in CCD coordinates).  p is the set of camera
% parameters: [dx,dy,psi0,th0,ph0,a].  The "badness of fit" vector
% is defined as the lengths of chords connecting catalog stars with stars
% in the image, when projected on a unit sphere.

[thccd,phccd] = cnv_xy2thph(star_x,star_y,p);
thcat = pi/2 - cat_el;
phcat = pi/2 - cat_az;
rccd = [sin(thccd).*cos(phccd);
        sin(thccd).*sin(phccd);
        cos(thccd)];
rcat = [sin(thcat).*cos(phcat);
        sin(thcat).*sin(phcat);
        cos(thcat)];
f = sqrt(sum((rccd-rcat).^2,1));
end

