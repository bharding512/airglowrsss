% Parameters
year = 2010;
fields = [5,10,15,20];
lon_bin = 5;
lat_bin = 5;

% Use AACGM to find magnetical lon/lat points in geographic coordinates and
% then trace along IGRF field lines to the conjugate hemisphere
la = 1:40;
lo = 0:359;
lat = la'*ones(size(lo));
lon = (lo'*ones(size(la)))';
[la,lo] = aacgm_cnv(lat,lon,zeros(size(lat)),'geo');

% Trace the IGRF to the conjugate hemisphere
for i = 1:length(la(:))
i
%     [lla,llaf] = IGRF_TraceFieldApex([la(i),lo(i),0],year,10e3);
%     [lla2,llaf2] = IGRF_TraceFieldConjugate(lla,year,0,10e3);
%     cla(i) = lla2(:,1);
%     clo(i) = lla2(:,2);
    % calculate the sign of the vertical component of the field
    [N,E,D,A] = igrfmx(0,year,1,0.1/1e3,90-la(i),lo(i));
    [lla,llaf] = IGRF_TraceField([la(i),lo(i),0.1],year,0,-sign(D).*10e3);
    cla(i) = lla(:,1);
    clo(i) = lla(:,2);

    % See if this is one of the field lines we want to save
    if((sum(lat(i) == fields) == 1) && (mod(lon(i),lon_bin) == 0))
        fieldlines(i) = {[llaf]};
    end
end

% Save the fieldlines
ge_plot_fieldlines(sprintf('MagneticFieldLines%d.kml',year),fieldlines);

% Loop through and save magnetic latitude grid
output = [];
for i = lat_bin:lat_bin:90
    outputx = ge_plot(lo(i,:),la(i,:),'LineWidth',2,'LineColor','ffffff00','name',sprintf('%2d N',i));
    output = [output, outputx];
    outputx = ge_plot(clo(i,:),cla(i,:),'LineWidth',2,'LineColor','ffffff00','name',sprintf('-%2d N',i));
    output = [output, outputx];
end

[mla,mlo] = aacgm_cnv(zeros(size(0:359)),0:359,zeros(size(0:359)),'geo');
outputx = ge_plot(lo(i,:),la(i,:),'LineWidth',2,'LineColor','ffff0000','name','Magnetic Equator');
output = [output,outputx];

ge_output('MagneticLatitude.kml',output);

% Loop through and save magnetic longitude grid
output = [];
for i = 1:lon_bin:359
    outputx = ge_plot([lo(:,i),clo(:,i)],[la(:,i),cla(:,i)],'LineWidth',2,'LineColor','ffffff00','name',sprintf('%2d E',i));
    output = [output, outputx];
end

ge_output('MagneticLongitude.kml',output);
