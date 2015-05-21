function [ccd_az, ccd_el, cost, p] = ...
        ImagerSpatialCal_BH(d, h, year, doy, time, p0, ... 
                 numAvailableStars, verbose)
%
% By Brian Harding 2012
%
% Match the stars in the image d to a known star field in order to infer
% the characteristics of the imager.  Use these characteristics to return
% the elevation and azimuth mappings. This is a more manual alternative to
% Ethan Miller's ImagerSpatialCal function.
%
% Inputs:
% d - MxN sky image
% h - image header, as returned by readtif (used for lat/lon of site)
% year - the year returned by tm2time(h.universaltime)
% doy - the day of year returned by tm2time(h.universaltime)
% time - the time of day (in sec) returned by
%       3600*tm2time(h.universaltime) PLUS the offset due to imperfect
%       timekeeping.  In short, this is when the beginning of exposure 
%       actually started, not just when it was reported to start.
% p0 - (1x7) The initial guess for the camera parameter vector (see below
%       for camera parameter vector definition).
% numAvailableStars - the number of stars in the catalog to use. (Typical:
%       100 for allsky, 1000 for narrow-field).
% verbose - If true, show the image, stars, and star field for visual
%           confirmation after convergence.
%
% Outputs:
% ccd_az - MxM mapping from pixel to azimuth, in degrees
% ccd_el - MxM mapping from pixel to elevation, in degrees
% cost - A scalar measure of how "far" the image stars are from the
%        known star field, in pixels.  Lower is better.
% p - The optimal imager parameter vector: [dx,dy,psi0,th0,ph0,a,flip]
%     (See below). This is probably not useful to the user.
%
% Definition of camera parameter vector: p = [dx,dy,psi0,th0,ph0,a,flip]
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
% a - rad/pix - the ratio between zenith angle and radial distance on the 
%   CCD. This parameter describes the mapping from real space to CCD
%   coordinates.
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
%   The first 6 parameters in the initial guess do not "need" to be close to
% the true values for the algorithm to converge, but they probably need to
% be close enough in order to manually select the correct stars.  The 7th
% parameter (the flip flag) needs to be correct for the algorithm to
% converge.
%   If the image is so noisy that star identification is impractical, try
% running starPassFilter to isolate stars before passing in an image.  Or
% subtract a dark frame.
%   If stars are still hard to differentiate from noise, try opening two
%   images and flipping back and forth.
% Example Code:
% fn = '../data/20101106_310/CAJ/6300_0180.tif';
% [d,header] = readtif(fn); % load image
% [year,~,~,doy,timehrs] = tm2time(header.universaltime);
% timesec = timehrs*3600; % convert to seconds
% % Parameters for Spatial Calibration
% numAvailableStars = 500; % Number of stars from the catalog to use
% verbose = 1; % If 1, show the result of convergence in an image.
% % These are just guesses based on the image.
% % initial guess of parameters: [dx,dy,del_psi,th0,ph0,a]
% dx = 250; % Center of image (x-coordinate)
% dy = 350; % Center of image (y-coordinate)
% psi0 = 0; % up is to the right ==> psi==0
% th0 = 80*pi/180; % 80 deg zenith angle
% ph0 = pi/2; % looking north ==> phi==pi/2
% a = (pi/6)/250; % It looks like 30 deg is about 250 pixels
% flipped = 1; 
% p0 = [dx,dy,psi0,th0,ph0,a,flipped];      
% % Run calibration
% [AZ,EL, ~, p] = ImagerSpatialCal_BH(d, header, year, doy, timesec,...
%         p0, numAvailableStars, verbose);
% % Save in the same format as airglowrsss/Matlab/SpatialCalibrations
% calibration = struct('az',AZ,'el',EL,'p',p);
% save('Calib','-struct','calibration');
%
% Future Work:
% 0) Allow the elevation falloff to be a polynomial (not just linear).
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


% (1) From star catalog, find elev/azim of the brightest stars.
    N = numAvailableStars; % Use the N brightest stars in the catalog
    el_cutoff = 0; % Stars in the catalog below this elevation angle will 
                    % not be used
    cat = load('yale_catalogue');
    [cat_az,cat_el] = cnv_cel2azel(year,doy,time,cat.obj_ra,...
        cat.obj_dec,h.latitude, h.longitude);
    m_el = find(cat_el > el_cutoff);
    [~,m_br] = sort(cat.obj_mag(m_el));
    m = m_el(m_br); % "m" indexes cat_az and cat_el
    cat_az = cat_az(m); % transpose to make column vector
    cat_el = cat_el(m);
    cat_br = cat.obj_mag(m);
    % Take N brightest
    cat_az = cat_az(1:N);
    cat_el = cat_el(1:N);
    cat_br = cat_br(1:N);
    % Convert to radians
    cat_az = (pi/180)*cat_az;
    cat_el = (pi/180)*cat_el;

% (2) Open a figure that has the image and the catalog stars
    h = figure;
    a1 = subplot(1,2,1);
    imagesc(d,[0 3*median(d(:))]);
    colormap('gray'); 
    axis equal
    axis tight;
    hold on;
    for r = [pi/6, pi/3, pi/2];
        tmp_ph = linspace(0,2*pi);
        tmp_th = 0*tmp_ph + r;
        [xxx,yyy] = cnv_thph2xy(tmp_th,tmp_ph,...
                p0(1),p0(2),p0(3),p0(4),p0(5),p0(6),p0(7));
        plot(xxx,yyy,'g-');
    end
    
    
    a2 = subplot(1,2,2);
%     imagesc(d,[0 3*median(d(:))]); hold on;
%     colormap('gray'); 
    cat_th = pi/2 - cat_el;
    cat_ph = pi/2 - cat_az;
    [xx,yy] = cnv_thph2xy(cat_th,cat_ph,...
                p0(1),p0(2),p0(3),p0(4),p0(5),p0(6),p0(7));
    s = 2.512.^(5.5 - cat_br);
    scatter(xx,yy,s,'k');
    set(gca,'YDir','Reverse');
    hold on;
    % plot circles
    for r = [pi/6, pi/3, pi/2];
        tmp_ph = linspace(0,2*pi);
        tmp_th = 0*tmp_ph + r;
        [xxx,yyy] = cnv_thph2xy(tmp_th,tmp_ph,...
                p0(1),p0(2),p0(3),p0(4),p0(5),p0(6),p0(7));
        plot(xxx,yyy,'k-');
    end
    hold off;
    axis equal
    axis tight
    xlim([0 size(d,2)]);
    ylim([0 size(d,1)]);
%     polar(cat_az,cat_ze_deg,'k.');
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
        %         % look around a little, so the user doesn't have to click right on
        %         % the star
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
    end
    cat_azi = cat_az(cat_idx);
    cat_eli = cat_el(cat_idx);

    close(h);
    pause(0.3); % to let it close
%     
    assignin('base','cat_azi',cat_azi);
    assignin('base','cat_eli',cat_eli);
    assignin('base','star_x',star_x);
    assignin('base','star_y',star_y);

% cat_azi = evalin('base','cat_azi');
% cat_eli = evalin('base','cat_eli');
% star_x = evalin('base','star_x');
% star_y = evalin('base','star_y');

% (4) Find the set of camera parameters p that optimizes star match.

%%%% gradient descent (slow) %%%%%
%     nm = 2;
%     f1 = @(pp) norm(star_fit_resid([pp p0(7)],...
%                         cat_azi,cat_eli,star_x,star_y),nm);
%                     
%     [p1,g] = gradientDescent(f1,p0(1:6),1e-8,1e-8);
%     assignin('base','g',g);
%     p = [p1(1:6) p0(7)];
%     cost = NaN;
                    
%%% fmincon %%%                    
    opts = optimset('fmincon');
    opts.Algorithm = 'interior-point'; % 'sqp' works fine too
    opts.MaxIter = 500;
    opts.TolX = 1e-6; % 1/12/12
    opts.TolFun = 1e-6;
    opts.Display = 'final';
%     opts.TolCon = 1e-8;
    
    % TODO: make this an input
    dp0 = [50 50 pi pi pi 0.01];
    lb = p0(1:6) - abs(dp0);
    ub = p0(1:6) + abs(dp0);
    nm = 2;
    f1 = @(pp) norm(star_fit_resid([pp p0(7)],...
        cat_azi,cat_eli,star_x,star_y),nm);
    
    [p1,cost,exitflag1,~,~,grad] = fmincon(f1,p0(1:6),[],[],[],[],lb,ub,[],opts);
    assignin('base','grad',grad);
    p = [p1(1:6) p0(7)];

%%% lsqnonlin %%%
%     scale = [100 100 1 1 1 0.01]; % scaling for numerical issues
%     p0(1:6) = p0(1:6)./scale;
%     f = @(pp) star_fit_resid([scale.*pp(1:6) p0(7)],...
%                             cat_azi,cat_eli,star_x,star_y);
%     opts = optimset('lsqnonlin');
%     opts.MaxIter = 50000;
%     opts.MaxFunEvals = 100000;
%     opts.TolX = 1e-6; % 1/12/12
%     opts.TolFun = 1e-6;
%     opts.Display = 'iter';
%     [p1,cost] = lsqnonlin(f,p0(1:6),[],[],opts);
%     p = [p1(1:6).*scale p0(7)];
    
% (5) Calculate and return the azimuth and elevation mappings from p.
    [X,Y] = meshgrid(1:size(d,2),1:size(d,1));
    [th,ph] = cnv_xy2thph(X(:)',Y(:)',p(1),p(2),p(3),p(4),p(5),p(6),p(7));
    ccd_el = zeros(size(d));
    ccd_az = zeros(size(d));
    ccd_el(:) = 90 - 180/pi * th;
    ccd_az(:) = 90 - 180/pi * ph;

% (6) If requested, plot the star fields
    if verbose
        figure;
        imagesc(d,[0 4*median(d(:))]); %scale image
        colormap gray
        hold on;
        plot(star_x,star_y,'go','MarkerSize',8);
        % TODO: plot more cat stars
        [cat_x,cat_y] = cnv_thph2xy(pi/2-cat_eli, pi/2-cat_azi, ...
                        p(1),p(2),p(3),p(4),p(5),p(6),p(7));
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
% TODO: flipped
dx = p(1);
dy = p(2);
psi0 = p(3);
th0 = p(4);
ph0 = p(5);
a = p(6);
flip = p(7);
[thccd,phccd] = cnv_xy2thph(star_x,star_y,dx,dy,psi0,th0,ph0,a,flip);
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

