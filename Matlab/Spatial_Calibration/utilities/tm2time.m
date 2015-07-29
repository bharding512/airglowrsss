function [year,month,day,doy,time] = tm2time(tm)
% TM2TIME converts the POSIX time struct 'tm' to human-readable form
%
% Syntax:
%
%   [year,month,day,doy,time] = TM2TIME(tm)
%
% Inputs:
%
%   tm - tm struct (from TIF header, for example)
%
% Outputs:
%   
%   year - year CE [years]
%   month - month [months, 1-12]
%   day - day of month [days, 1-31]
%   doy - day of year [days, 1-366]
%   time - time [hrs]
%
% TM2TIME was written by:
%
% Ethan Miller (esmiller@illinois.edu)
% Remote Sensing and Space Sciences Group (http://rsss.csl.uiuc.edu/)
% Department of Electrical and Computer Engineering
% University of Illinois at Urbana-Champaign
% Current:  2009 April 16 / Version 1.0
% Original:  2009 April 16 / Version 1.0     

    doy = 1+tm(15);
    year = 1900+tm(11);
    month = 1+tm(9);
    day = tm(7);
    time = tm(5)+tm(3)/60+tm(1)/3600;
    
    