function MST = UTC2MST(lon,doy,year,utc)
% UTC2MST computes the Mean Sidereal Time (MST) for a given UTC and
% longitude
%
% Syntax:
%
%   MST = UTC2MST(lon,doy,year,utc)
%
% Inputs:
%
%   lon - east longitude [deg]
%   doy - day-of-year [days]
%   year - year [years]
%   utc - universal time [sec]
%
% Outputs:
%   
%   MST - Mean Sidereal Time [deg]
%
% UTC2MST was written by:
%
% Ethan Miller (esmiller@illinois.edu) and Peter Hedlund (hedlundp@gmail.com)
% Remote Sensing and Space Sciences Group (http://rsss.csl.uiuc.edu/)
% Department of Electrical and Computer Engineering
% University of Illinois at Urbana-Champaign
% Current:  2009 April 16 / Version 1.0
% Original:  2009 April 16 / Version 1.0 

    % Calculate the sidereal constants needed to convert from UTC to GMST.
    % This is from http://celestrak.com/columns/v02n01/
    dd = doy;
    yy = year;
    hh = floor(utc/3600);
    mm = floor((utc/3600-hh)*60.);
    ss = round(((utc/3600-hh)*60.-mm)*60.);
    jd = date2jd(yy,1,dd,hh,mm,ss)-2451545;
    
    % This formula is from:
    % http://www2.arnes.si/~gljsentvid10/sidereal.htm
    % Greenwich MST
    GMST = 280.46061837 + 360.98564736629 * jd;
    % Local MST
    MST = mod(GMST+lon,360);

% Alternate code that I believe gets the same result from:
% http://home.att.net/~srschmitt/celestial2horizon.html
%
%     dm = doy2date(doy,year);
%     day = dm(:,1);
%     month = dm(:,2);
%     
%     a = floor( year/100 );
%     b = 2 - a + floor( a/4 );
%     c = floor( 365.25*year );
%     d = floor( 30.6001*( month + 1 ) );
% 
%     % days since J2000.0
%     jd = b + c + d - 730550.5 + day + (hh + mm/60.0 + ss/3600.0)/24.0;
%     
%     % julian centuries since J2000.0
%     jt = jd/36525.0;
% 
%     % mean sidereal time
%     MST = 280.46061837 + 360.98564736629*jd ...
%          + 0.000387933*jt*jt - jt*jt*jt/38710000 + lon;

