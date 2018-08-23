function [trise,tset,daylen]=shinetime(lon0,lon,lat,yy,mm,dd)
clear global; clear functions;
persistent cost delta eta jday jday0 pi t w ; 
%#############################################################################;
%     purpose:;
%     calculate maximum possible sun shine time(hour);
%#############################################################################;

%input:;
% the longitude where the standard time is defined.
% If input time is UTC, lon0 = 0
% If input time is BST, lon0 = 120
% If input time is JST, lon0 = 135
% longitude of each site (deg.)
% > 0 in East Hemesphere
% < 0 in West Hemesphere
% latitdue of each site (deg.)
% > 0 in North Hemesphere
% < 0 in South Hemesphere
%calendar year:month:day

%output:;
% sunrise (hr)
% sunset  (hr)
% sunset-sunrise
if isempty(pi), pi=3.1415926; end;
% day number starting with 1 on Jan.1
if isempty(jday), jday=0; end;
if isempty(jday0), jday0=0; end;
% solar declination(rad)
if isempty(delta), delta=0; end;
if isempty(eta), eta=0; end;
if isempty(t), t=0; end;
if isempty(cost), cost=0; end;
if isempty(w), w=0; end;
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@;
%;
%Begining of executable code......
%;
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@;
[jday] = YYYYMMDD2doy(yy, mm, dd);
[jday0] = YYYYMMDD2doy(yy, 12, 31);

w    = fix(((2.*pi).*jday)./jday0);

%calculate solar declination angle delta from calendar;
delta= 0.3622133-23.24763*cos(w+0.153231)-0.3368908*cos(2*w+0.2070988)-0.1852646*cos(3*w+0.6201293);
delta= delta * pi/ 180;

%calculate daily averaged time error due to elliptic size of earth;
eta  = 60*(-0.0002786409+0.1227715*cos(w+1.498311)- 0.1654575*cos(2*w-1.261546)-0.00535383*cos(3*w-1.1571));
cost=-tan(delta).*tan(lat*pi/180);
cost=max(-1,min(1,cost));
t=acos(cost)/pi*180;

% (hour)
trise=max(0,(-t-(lon-lon0))/15-eta./60+12);
% (hour)
tset=min(24,( t-(lon-lon0))/15-eta./60+12);

daylen = tset-trise;

return;
end %subroutine shinetime

