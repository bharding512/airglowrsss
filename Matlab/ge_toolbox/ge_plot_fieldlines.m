function ge_plot_fieldlines(gefilename,fieldlines, date_start, date_stop)
% The GE_PLOT_FIELDLINES function plots geomagnetic field lines to
% Google Earth KML format.
%
% Syntax:
%
%   GE_PLOT_FIELDLINES(FILENAME,FIELDLINES)
%
% Inputs:
%
%   FILENAME - output filename string
%   FIELDLINES - cell array field line Lat, Lon, Alt [deg,deg,m] matrices
%
% Outputs:
%
%   none
%
% TODO:
%
%   A newer version will allow the line style to be passed.
%
% GE_PLOT_FIELDLINES was written by:
%
% Ethan Miller (esmiller@illinois.edu)
% Remote Sensing and Space Sciences Group (http://rsss.csl.illinois.edu/)
% Department of Electrical and Computer Engineering
% University of Illinois at Urbana-Champaign
% Current:  2009 February 5 / Version 1.0 (original part of COSMICscintsGE)

output = [];
for i=1:length(fieldlines)
    lla = cell2mat(fieldlines(i));
    if(length(lla) > 1)
        clear outputx;
        outputx = ge_plot3(lla(:,2),lla(:,1),lla(:,3),'LineWidth',2,...
            'LineColor','ffffff00','timeSpanStart',date_start,'timeSpanStop',date_stop);
        output = [output outputx];
    end
end

ge_output(gefilename,output);
