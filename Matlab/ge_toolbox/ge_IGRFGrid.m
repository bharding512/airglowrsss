function ge_IGRFGrid(year, fields_mlat, delta_mlon)
%
% A function that generates KML files containing IGRF approximations of the
% magnetic field.  The function will generate magnetic field lines that
% begin from the magnetic latitudes requested.  These locations are
% calcualted from IGRF by finding the angle between the
% IGRF main field at the ground and the line from the center of the Earth
% and the point on the ground.  This is done, rather than using AACGM, as
% AACGM is no longer maintained.
%
% Inputs:
%   year - year to run IGRF
%   fields_mlat - a vector containing the magnetic latitudes (N) for which
%                 field lines should be calculated
%   delta_mlon - the longitudinal spacing between field lines
%
% Written by Jonathan J. Makela on 4 Feb 2010
%

% Define vectors
geo_lon = 0:1:360;
geo_lat = -40:1:40;

% String for the time stamp
date_start = sprintf('%04d-01-01T00:00', year);
date_stop = sprintf('%04d-12-31T23:59', year);

% Calculate lines of constant magnetic latitude
output_mlat = [];
nl = zeros([length(geo_lon),3]);
sl = zeros([length(geo_lon),3]);
for lat = fields_mlat
    for j = 1:length(geo_lon)
        disp(sprintf('1: %s Lat: %5.1f Lon: %5.1f', datestr(now), lat, geo_lon(j)))
        [nlla,slla] = igrf_mag2geo(lat,geo_lon(j),0,year);
        nl(j,:) = nlla;
        sl(j,:) = slla;
    end
    if(lat == 0)
        output = ge_plot(nl(:,2),nl(:,1),'LineWidth',2,'LineColor','ffffff00','name','Magnetic Equator','timeSpanStart',date_start,'timeSpanStop',date_stop);
        output_mlat = [output_mlat, output];
    else
        output = ge_plot(nl(:,2),nl(:,1),'LineWidth',2,'LineColor','ffffff00','name',sprintf('%2d N',lat),'timeSpanStart',date_start,'timeSpanStop',date_stop);
        output_mlat = [output_mlat, output];
        output = ge_plot(sl(:,2),sl(:,1),'LineWidth',2,'LineColor','ffffff00','name',sprintf('%2d S',lat),'timeSpanStart',date_start,'timeSpanStop',date_stop);
        output_mlat = [output_mlat, output];
    end
end

% Write the file
ge_output(sprintf('MagneticLatitude_%d.kml', year), output_mlat);

%Calculate lines of constant magnetic longitude
output_mlon = [];
nl = zeros([length(geo_lat),3]);
sl = zeros([length(geo_lat),3]);
for lon = 0:delta_mlon:360
    for j = 1:(length(geo_lat)/2+1)
        disp(sprintf('2: %s Lat: %5.1f Lon: %5.1f', datestr(now), geo_lat(j), lon))
        [nlla,slla] = igrf_mag2geo(geo_lat(j),lon,0,year);
        ind2 = find(nlla(:,2) > 180);
        nlla(ind2,2) = nlla(ind2,2) - 360;
        ind2 = find(slla(:,2) > 180);
        slla(ind2,2) = slla(ind2,2) - 360;
        nl(j,:) = nlla;
        nl(length(sl)+1-j,:) = slla;
    end
    
    output = ge_plot(nl(:,2),nl(:,1),'LineWidth',2,'LineColor','ffffff00','name',sprintf('%2d E',lon),'timeSpanStart',date_start,'timeSpanStop',date_stop);
    output_mlon = [output_mlon, output];
end

% Write the file
ge_output(sprintf('MagneticLongitude_%d.kml', year), output_mlon);

% Calculate field lines
lines_count = 0;
output_fields = [];
for lat = fields_mlat
    for lon = 0:delta_mlon:360
        disp(sprintf('3: %s Lat: %5.1f Lon: %5.1f', datestr(now), lat, lon))
        [nlla,slla] = igrf_mag2geo(lat,lon,0.1,year);
        nlla(3) = 0.1;
        slla(3) = 0.1;
        [E,N,D,F] = igrf(year,nlla(1),nlla(2),nlla(3));
        [lla,llaf] = IGRF_TraceField(nlla,year,0,-sign(D).*10e3);
        
        ind2 = find(llaf(:,2) > 180);
        llaf(ind2,2) = llaf(ind2,2) - 360;
        lines_count = lines_count+1;
        fieldlines(lines_count) = {[llaf]};
    end
end

output_fields = [output_fields, fieldlines];
ge_plot_fieldlines(sprintf('MagneticFieldLines_%d.kml', year), output_fields, date_start, date_stop);