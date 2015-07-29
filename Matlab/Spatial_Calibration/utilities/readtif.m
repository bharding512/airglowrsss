% [data, header] = readtif(fname)
%
% Function to read in a TIF file and return the data as well as the header
% file.  Note, this should correctly read in both old (GMS) TIF files and
% new (CDF) TIF files.
%
% INPUTS:
%   fname - path and filename of file to be read.
%
% OUTPUTS:
%   data - the intensity data saved in the file.
%   header - the header information saved in the file.
%
% HISTORY:
%   Sep 23, 2005 - Written by Jonathan J. Makela based on historic code
%   written for IDL and C++

function [data, header] = readtif(fname);

% Check to make sure the image exists
if(exist(fname) ~= 0)
    % open up the file for reading
    fid = fopen(fname, 'r');
    
    % skip the TIF header stuff
    fseek(fid,200,'bof');
    
    % read in the image type (should be 'CDF')
    header.type = sprintf('%s',char(fread(fid,3,'uchar')));
    if(strcmp(header.type,'CDF'))
        isGMS = 0;
    else
        isGMS = 1;
    end
    
    % read in the number of rows and columns
    header.rows = fread(fid,1,'uint16','ieee-le');
    header.cols = fread(fid,1,'uint16','ieee-le');
    
    % read in the autoscaled min and max pixel values
    header.min = fread(fid,1,'uint16','ieee-le');
    header.max = fread(fid,1,'uint16','ieee-le');
    
    % read in the local and universal time structures
    header.localtime = fread(fid,18,'uint16','ieee-le');
    header.universaltime = fread(fid,18,'uint16','ieee-le');
    
    % read in the filter name
    header.filter = sprintf('%s',char(fread(fid,12,'uchar')));
    
    % read in the emission height
    header.emissionheight = fread(fid,1,'float','ieee-le');
    
    % read in the exposure time
    header.exposuretime = fread(fid,1,'float','ieee-le');
    
    % read in CCD parameters
    header.gain = fread(fid,1,'uchar');
    header.XBinning = fread(fid,1,'uchar');
    header.YBinning = fread(fid,1,'uchar');
    header.CCDTemperature = fread(fid,1,'float','ieee-le');
    
    % read in FW temperature
    header.FWTemperature = fread(fid,1,'float','ieee-le');
    
    % read in XY spacing from unwarping
    header.XSpacing = fread(fid,1,'float','ieee-le');
    header.YSpacing = fread(fid,1,'float','ieee-le');
    
    % read in the spatial calibration coefficients
    header.calibration.a0 = fread(fid,1,'double','ieee-le');
    header.calibration.a1 = fread(fid,1,'double','ieee-le');
    header.calibration.a2 = fread(fid,1,'double','ieee-le');
    header.calibration.b0 = fread(fid,1,'double','ieee-le');
    header.calibration.b1 = fread(fid,1,'double','ieee-le');
    header.calibration.b2 = fread(fid,1,'double','ieee-le');
    
    % read in the site location
    header.latitude = fread(fid,1,'float','ieee-le');
    header.longitude = fread(fid,1,'float','ieee-le');
    header.altitude = fread(fid,1,'float','ieee-le');
    header.observatory = sprintf('%s',char(fread(fid,30,'uchar')));
    
    % old images are referenced to west longitude, not east
    if(isGMS)
        header.longitude = 360 - header.longitude;
    end
    
    % read in the comment
    header.comment = sprintf('%s',char(fread(fid,100,'uchar')));
    
    % the CDF files have some additional information
    if(~isGMS)
        header.isAllSky = fread(fid,1,'uchar');
        header.CenterAz = fread(fid,1,'double','ieee-le');
        header.CenterEl = fread(fid,1,'double','ieee-le');
        header.ProjectionAltitude = fread(fid,1,'float','ieee-le');
        header.ProjectionLon = fread(fid,1,'float','ieee-le');
        temp = fread(fid,63,'uchar');
    else
        temp = fread(fid,88,'uchar');
    end
    
    % read in the data
    data = fread(fid,header.rows*header.cols,'ushort','ieee-le');
    data = reshape(data,header.cols,header.rows);
    
    % close the file
    fclose(fid);
end
    
    