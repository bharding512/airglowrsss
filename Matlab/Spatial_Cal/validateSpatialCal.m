function validateSpatialCal( im_fn, cal_fn )
% Display an image to visually validate the saved spatial calibration for a
% desired image.
% Inputs:
%   im_fn: The filename of the desired image (*.tif only)
%   cal_fn: The filename of the spatial calibration file
% Output:
%   Pop up a figure showing the image with calibrated star positions shown
%   with green circles.  If they line up with the image stars, then the
%   spatial calibration is validated.
%
% Example Call:
%  im_fn = ...
%  '/home/harding/airglowrsss/Harding/Imager_Tomography/data/20101106_310/CAJ/6300_0075.tif';
%  cal_fn = ...
%  '/home/harding/airglowrsss/Matlab/ImagerCalibrations/PICASSO5elaz_2652009';
%  validateSpatialCal(im_fn, cal_fn);

el_thresh = 20; % Threshold above which to show stars
N = 50; % Number of stars to show

% load stuff
    [im,header] = readtif(im_fn); % image
    calib = load(cal_fn); % calibration
    EL = calib.el;
    AZ = calib.az;

% obtain current time
    [yr,~,~,doy,ut] = tm2time(header.universaltime);
  
% get catalogue star locations    
    cat = load('yale_catalogue');
    [cat_az,cat_el] = cnv_cel2azel(yr, doy, ut*3600, cat.obj_ra,...
                cat.obj_dec, header.latitude, header.longitude);
    m_el = find(cat_el > el_thresh); % Only use stars above a threshold
    [~,m_br] = sort(cat.obj_mag(m_el));
    m = m_el(m_br); % "m" indexes cat_az and cat_el
    cat_az = cat_az(m);
    cat_el = cat_el(m);

% For each star, find the pixel closest to the star position
    stari = zeros(N,1);
    starj = zeros(N,1);
    for i = 1:N
        e = cat_el(i);
        a = cat_az(i);
        [~,idx] = min(abs(EL(:) - e) + abs(AZ(:) - a));
        [si,sj] = ind2sub(size(EL),idx);
        stari(i) = si;
        starj(i) = sj;
    end

% Plot
    imagesc(im,[0 prctile(im(:),99.9)]);
    colormap('gray');
    hold on;
    plot(starj,stari,'go','MarkerSize',5);
    hold off;


end

