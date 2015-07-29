function [az,el] = cnv_cel2azel(year,doy,utc,RA,dec,lat,lon)
% CNV_CEL2AZEL computes the azimuth and elevation of an object given its 
% right-ascension and declination and an observing location/date/time
%
% Syntax:
%
%   [az,el] = cnv_cel2azel(year,doy,utc,RA,dec,lat,lon)
%
% Inputs:
%
%   year - year [years]
%   doy - day-of-year [days]
%   utc - universal time [sec]
%   RA - right ascension [deg]
%   dec - declination [deg]
%   lat - observing north latitude [deg]
%   lon - observing east longitude [deg]
%
% Outputs:
%   
%   az - azimuth [deg]
%   el - elevation [deg]
%
% CNV_CEL2AZEL was written by:
%
% Ethan Miller (esmiller@illinois.edu) and Peter Hedlund (hedlundp@gmail.com)
% Remote Sensing and Space Sciences Group (http://rsss.csl.uiuc.edu/)
% Department of Electrical and Computer Engineering
% University of Illinois at Urbana-Champaign
% Current:  2009 April 16 / Version 1.0
% Original:  2009 April 16 / Version 1.0 
%
% based on code from:
%
% http://home.att.net/~srschmitt/celestial2horizon.html

    % compute hour angle in degrees
    HA = UTC2MST(lon,doy,year,utc) - RA;
    ix = find(HA < 0);
    HA(ix) = HA(ix) + 360;
    
    % convert degrees to radians
    HA  = HA*pi/180;
    dec = dec*pi/180;
    lat = lat*pi/180;

    % compute elevation in radians
    sin_alt = sin(dec).*sin(lat) + cos(dec).*cos(lat).*cos(HA);
    el = asin(sin_alt);
    
    % compute azimuth in radians
    % divide by zero error at poles or if alt = 90 deg
    cos_az = (sin(dec) - sin(el).*sin(lat))./(cos(el).*cos(lat));
    az = acos(cos_az);

    % convert radians to degrees
    el = el*180/pi;
    az = az*180/pi;

    % choose hemisphere
    ix = find(sin(HA) > 0);
    az(ix) = 360 - az(ix);
