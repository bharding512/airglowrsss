function monstr = month2str(month,format)
% The MONTH2STR function converts an enumerated month to its name. 
%
% Syntax:
%
%   MONSTR = MONTH2STR(MONTH,[FORMAT])
%
% Inputs:
%
%   MONTH - Month, enumerated in [1,12]
%   FORMAT - optional format, if set to 'long', long names will be used
%            any other input (including excluding it) returns abreviations
%
% Outputs:
%   
%   MONSTR - cell array of month names
%
% MONTH2STR was written by:
%
% Ethan Miller (esmiller@uiuc.edu)
% Remote Sensing and Space Sciences Group (http://rsss.csl.uiuc.edu/)
% Department of Electrical and Computer Engineering
% University of Illinois at Urbana-Champaign
% Current:  2009 February 3 / Version 1.0
% Original:  2009 February 3 / Version 1.0  

    monthlong = {'January','February','March','April','May','June',...
        'July','August','September','October','November','December'};
    monthshort = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep',...
        'Oct','Nov','Dec'};
    
    for i = 1:length(month)
        if ((nargin==2) & (format == 'long'))
            monstr(i) = monthlong(month(i));
        else
            monstr(i) = monthshort(month(i));
        end
    end