function [loz]=NASA_ozone(year,jday,lat,lon)
%#######################################################################;
%     purpose:;
%     calculate the ozone thickness according to the monthly dataset of;
%     nasa toms;
%#######################################################################;
%  input:
% jday: julian day
% lat: latitude
% lon: longitude
% Local latitude, positive for northern hemisphere (deg.)
% longitude of each site (deg.)
% > 0 in East Hemesphere
% < 0 in West Hemesphere

%  output: ozone thickness

load ozone_monthly_nasa;
[month day] = doy2MMDD(year, jday);
rr=(lat>=ozone_monthly_nasa(:,1))&(lat<=ozone_monthly_nasa(:,2));
loz = ozone_monthly_nasa(rr,month+2);
return;

