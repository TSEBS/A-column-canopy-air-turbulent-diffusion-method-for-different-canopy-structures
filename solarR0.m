function [R0,hsun]=solarR0(lon0,lon,lat,yy,mm,dd,hh,mn,ss)
%#############################################################################;
%
%     purpose: calculate horizontal extraterrestial solar insolation(w./s);
%
%#############################################################################;
% input:;
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
% calendar year:month:day
% standard time hour:minute:second in lon0

% output:;
% R0: horizontal extraterrestial solar insolation(w/m^2) 
% hsun: the height of the sun (rad)
% solar constant
i00=1367; 
pi=3.1415926;
%(a^2/r^2):earth-sun distance factor

%Begining of executable code......
[jday] = YYYYMMDD2doy(yy, mm, dd);
[jday0] = YYYYMMDD2doy( yy, 12, 31);

w    = (2*pi*jday/jday0); % Kumar 1997, equation 7
%     calculate square distance between the earth and the sun:(d0./d)^2;
d0d2 = 1.00011+0.034221*cos(w)+0.00128*sin(w)+ 0.000719*cos(2*w)+0.000077*sin(2*w);% Bob`s method %  Kumar 1997, equation 7

%calculate solar declination angle delta from calendar;
delta= 0.3622133-23.24763.*cos(w+0.153231)-0.3368908.*cos(2*w+0.2070988)-0.1852646.*cos(3*w+0.6201293);
delta= delta * pi / 180;

%calculate daily averaged time error due to elliptic size of earth orbit(min);
eta  = 60*(-0.0002786409+0.1227715*cos(w+1.498311)-0.1654575*cos(2*w-1.261546)-0.00535383 * cos(3*w-1.1571));



ts   = (hh)+(mn)./60+(ss)./3600; %convert standard time
%     convert hour and longitude difference to hour angle(rad);
t    = 15.*(ts-12 +eta/60)+lon-lon0;
hrangle=t* pi/180; % hour angle(rad)
%calculate the height of the sun(rad);
sinh  = sin(lat.*pi./180).*sin(delta) + cos(lat.*pi./180).*cos(delta).*cos(hrangle);

% quality control
hsun=sinh;
hsun(sinh<=0)=0;
hsun(sinh>0)=asin(sinh(sinh>0));
R0 = i00 .* d0d2 .* sin(hsun);

return;
end

