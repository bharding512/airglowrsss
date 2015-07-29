%% Example code for using ImagerSpatialCalib

% Load allsky image and other information
    fn = '../066/6300_0015.tif';
    d = readtif(fn); % the matrix of intensity values (the image)
    lat = -6.8710; % latitude of site
    lon = 321.4420; % longitude of site
    year = 2010; % year the image was taken
    doy = 311; % day of year the image was taken (1-366)
    timesec = 15602; % time of day in seconds

% Initial approximate guess for camera parameters
    dx = size(d,2)/2; % center pixel x-coordinate
    dy = size(d,1)/2; % center pixel y-coordinate
    psi0 = 0; % for allsky images, controls the rotation of the image (rad)
    th0 = 0; % for allsky images, set this to zero
    ph0 = 0; % for allsky images, set this to zero
    a = [pi/size(d,2),0,0]; % coefficients of cubic function
    flipped = 0; % this controls front- vs. back-illuminated CCD (0 or 1)
    p0 = [dx,dy,psi0,th0,ph0,a,flipped]';
%     p0 = p; % use the last solved vector as the initial guess
    % a priori known flag
    known = [0 0 0 0 0 0 0 0 0]';

% Run Calibration
    numAvailableStars = 200; % controls how many catalog stars to display
    verbose = 1; % if 1, displays output image
    useold = 0;
    [ccd_az, ccd_el, cost, p] = ImagerSpatialCalib(d, lat, lon, year, ...
                    doy, timesec, p0, numAvailableStars, verbose,...
                    useold, known);


%% Use the following sections if you have a directory of data with the
%  header structure used in Makela's group.
%% Load Images from a directory
% Create variables to use later: 
%       fns, years, doys, times, sitelatlon, numImages

time_offset = 0;

folder = '../picasso02_cto_2015209';
fns = cell([1000,1]);
times = [];
years = [];
doys = [];

d = dir(folder);
i=1;
for k=1:length(d)
    if ~isempty(strfind(d(k).name,'.tif')) && ~isempty(strfind(d(k).name, '6300'))
        fn = [folder '/' d(k).name];
        [~,header] = readtif(fn);
        [yr,mo,day,doy,ut] = tm2time(header.universaltime);
        utsec = ut*3600;
        fns{i} = fn;
        times(i) = mod(utsec - time_offset,86400);
        % figure out doy
        if (utsec - time_offset > 86400)
            realdoy = doy + 1; % TODO: year rollover
        elseif (utsec - time_offset < 0)
            realdoy = doy - 1;
        else
            realdoy = doy;
        end
        doys(i) = realdoy;
        years(i) = yr;
        
        i=i+1;
    end
end
sitelatlon(1) = header.latitude;
sitelatlon(2) = header.longitude;

numImages = i-1;

%% Run Spatial Cal
% (Don't forget to change the filename of the saved calibration file below)
% (Hint: on Mac, use Cmd+tilde to switch between windows)

% Optional: specify dark image
    usedark = 1;
    darkfn = '../picasso02_cto_2015209/CTO_DARK_20150728_093127_0003.tif';

% Load allsky image and other information
    i = 33;
    fn = fns{i};
    d = readtif(fn); % the matrix of intensity values (the image)
    if usedark
        dk = readtif(darkfn);
        d = d - dk;
    end
    sitelat = sitelatlon(1);
    sitelon = sitelatlon(2);
    year = years(i);
    doy = doys(i);
    timesec = times(i);

% Initial approximate guess for camera parameters
    dx = size(d,2)/2; % center pixel x-coordinate
    dy = size(d,1)/2; % center pixel y-coordinate
    psi0 = 0.0; % for allsky images, controls the rotation of the image (rad)
    th0 = 0.95*pi/2; % for allsky images, set this to zero (ze of center pixel)
    ph0 = pi/2; % for allsky images, set this to zero (az of center pixel)
    a = [0.5*pi/size(d,2),0,0]; % coefficients of cubic function
    flipped = 1; % this controls front- vs. back-illuminated CCD (0 or 1)
    p0 = [dx,dy,psi0,th0,ph0,a,flipped]';
    p0 = p; % use the last solved vector as the initial guess
    % a priori known flag
    known = [0 0 0 0 0 0 0 0 0]';
    
% Run Calibration
    numAvailableStars = 500; % controls how many catalog stars to display
    verbose = 1; % if 1, displays output image
    useold = 0;
    [ccd_az, ccd_el, cost, p] = ImagerSpatialCalib(d, sitelat, sitelon, year, ...
                    doy, timesec, p0, numAvailableStars, verbose,...
                    useold, known);

% Save Calibration
    % generate lat/lon projected maps
    H = 250;
    [lat lon] = cnv_azel2latlon(ccd_az,ccd_el,H,sitelatlon);
    calibration = struct('az',ccd_az,'el',ccd_el,'lat',lat,...
                         'lon',lon,'p',p,'sitelat',sitelat,'sitelon',...
                         sitelon);
    save('CTOelaz_2015209','-struct','calibration');
    
    
%% Visually confirm the saved calibration by showing movie

% Optional: specify dark image
    usedark = 1;
    darkfn = '../picasso02_cto_2015118/CTO_DARK_20150428_075944_0006.tif';

% Load spatial cal
    calib = load('CTOelaz_2015209.mat');
    
% folder to save images
    folder = sprintf('mor_figures/');

% now make images
    close all;
    figure;
    set(gcf,'Position',[1357   456   1217    1019]);
    
    for i = 1:numImages
        % load data
        [d,h] = readtif(fns{i});
        if usedark
            dk = readtif(darkfn);
            d = d - dk;
        end

        % get catalogue star locations
        N = 500;
        cat = load('yale_catalogue');
        year = years(i);
        doy = doys(i);
        time = times(i);
        [cat_az,cat_el] = cnv_cel2azel(year,doy,time,cat.obj_ra,...
                    cat.obj_dec,h.latitude, h.longitude);
        m_el = find(cat_el > 5);
        [~,m_br] = sort(cat.obj_mag(m_el));
        m = m_el(m_br); % "m" indexes cat_az and cat_el
        cat_az = cat_az(m); % transpose to make column vector
        cat_el = cat_el(m);   
        cat_az = cat_az(1:N);
        cat_el = cat_el(1:N);
        
        % Convert to ccd coordinates
        [x,y] = cnv_thph2xy(pi/180*(90-cat_el),pi/180*(90-cat_az),calib.p);
        
        % Plot
        imagesc(d,prctile(d(:),[40,99.9]));
        colormap('gray');
        hold on;
        plot(x,y,'m^','MarkerSize',8);
        hold off;
        title(i)
        pause(0.1);
        
        
        % Save image
%         export_fig([folder sprintf('%04i',i)],'png');
        %saveas(gcf,[folder sprintf('%04i',i) '.png']);
    end
    
%% Show Movie of Data (so that you can pick the best image)

usedark = 1;
darkfn = '../picasso02_cto_2015118/CTO_DARK_20150428_075944_0006.tif';

for i = 1:numImages
    d = readtif(fns{i});
    if usedark
        dk = readtif(darkfn);
        d = d - dk;
    end
    imagesc(d,prctile(d(:),[15 99])); colormap gray;
    title(i)
    pause(0.01);
end
    
    
    
    
