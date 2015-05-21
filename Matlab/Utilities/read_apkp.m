function [dn, kp, daily_kp, ap, mean_ap, f107] = read_apkp(year)
%
% Reads in the apkp formatted data as described in kp_ap.fmt.  Data
% downloaded from ftp.ngdc.noaa.gov/STP/GEOMAGNETIC_DATA/INDICES/KP_AP/
%
% INPUT
%   year - the year to read in
%
% OUTPUT
%   dn - a [1,nday] array containing the datenum for each day
%   kp - a [8,nday] array containing the 3-hour Kp index
%   daily_kp - a [1,nday] array containing the sum of the Kp for the day
%   ap - a [8,nday] array containing the 3-hour ap index
%   mean_ap - a [1,nday] array containing the mean ap for the day
%   f107 - a [1,nday] array containing the Ottawa 10.7-cm solar radio flux
%
% Written by Jonathan J. Makela on 12 Feb 2008 (jmakela@uiuc.edu)
% Modified by Jonathan J. Makela on 29 Mar 2010 to calculate the path_stub

% Find the location of this function and then create the path to the data
% files based upon it.
temp = which('read_apkp.m');
path_stub = sprintf('%s/ap_values/', temp(1:end-11));

% Open the requested file
fname = sprintf('%s%d', path_stub, year);
fid = fopen(fname);

% Read in the data
%d = fscanf(fid,'%2d%2d%2d%4d %2d %2d %2d %2d %2d %2d %2d %2d %2d %3d %3d %3d %3d %3d %3d %3d %3d %3d %5f %1d %3d %6f',[26,inf]);
i = 1;
valid = 1;
line = fgetl(fid);
while valid
    d(:,i) = [str2double(line(1:2)), str2double(line(3:4)),...
        str2double(line(5:6)), str2double(line(7:10)),...
        str2double(line(11:12)), str2double(line(13:14)),...
        str2double(line(15:16)), str2double(line(17:18)),...
        str2double(line(19:20)), str2double(line(21:22)),...
        str2double(line(23:24)), str2double(line(25:26)),...
        str2double(line(27:28)), str2double(line(29:31)),...
        str2double(line(32:34)), str2double(line(35:37)),...
        str2double(line(38:40)), str2double(line(41:43)),...
        str2double(line(44:46)), str2double(line(47:49)),...
        str2double(line(50:52)), str2double(line(53:55)),...
        str2double(line(56:58)), str2double(line(62)),...
        str2double(line(63:65)), str2double(line(66:70))]';

    line = fgetl(fid);
    i = i+1;
    if(line == -1)
        valid = 0;
    end
end


% Close the file
fclose(fid);

% Parse the data
y = year * ones(1,length(d));
m = d(2,:);
da = d(3,:);
dn = datenum(y,m,da);
kp = d(6:13,:);
daily_kp = d(14,:);
ap = d(15:22,:);
mean_ap = mean(ap,1);    %d(23,:);
f107 = d(26,:);

end