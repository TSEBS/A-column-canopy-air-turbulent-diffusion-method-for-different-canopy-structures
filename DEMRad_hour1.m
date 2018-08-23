function [R0,Rsfc,Rb,Rd,Rr,taub,taud,taur]=DEMRad_hour1(lon,lat,yy,mm,dd,hh,mn,ss,pa,ta,rh,DEM,albedo,sunrs_switch,shine)
%  INPUT:
%      real    lon        ! longitude of each site (deg.)
%                         ! > 0 in East Hemesphere
%                         ! < 0 in West Hemesphere
%      real    lat        ! latitdue of each site (deg.)
%                         ! > 0 in North Hemesphere
%                         ! < 0 in South Hemesphere
%      integer	yy,mm,dd,hh,hh, mn,ss, GMT time
%      real    pa         ! surface air pressure (Pa)
%      real    ta         ! air temperature (K)
%      real    rh         ! relative humidity
%      real    shine      ! Hourly sunshine duration (hour)
%      real    DEM

%  OUTPUT:
%     real    Rsfc        ! Hourly all-sky solar irradiance (W m^-2)

%  Temporary:
%      integer jday       ! julian day
%      integer	mn,ss     ! standard time minute:second in lon0
%      real    R0         ! Instantaneous TOA solar irradiance (W m^-2)
%      real    radb       ! Hourly clear-sky solar beam (W m^-2)
%      real    radd       ! Hourly clear-sky solar diffuse (W m^-2)
%      real    rad        ! Hourly clear-sky solar irradiance (W m^-2)
%      real    shinemax   ! Hourly maximum sunshine duration (hour)
%      real hsun      ! the height of the sun (90-zenith angle)
%      real taub      ! Beam radiative transmittance
%      real taud      ! Diffuse radiative transmittance
%      real ratio     ! realtive sunshine duration 
%      real trise     ! sunrise time (hour)
%      real tset      ! sunset  time (hour)
%      real daylen    ! sunset-sunrise
%      real dt        ! time step for integrating hourly solar radiation (hour)

lon0   = 0;
[jday] = YYYYMMDD2doy(yy, mm, dd);
tau_a  = yeardays(yy);     % length of the year in days
%     calculate sunrise and sunset in hours
[m,n] = size(DEM);
dt   = [yy mm dd];
% rs    = zeros(size(DEM)); 

if sunrs_switch==1 % the following code to compute sunrise and sunset, when used to process global satellite dataset
    for i=1:m
        for j=1:n
            [rs(i,j,:)] = suncycle(lat(i,j),lon(i,j),dt);   % output rs in GMT time
        end
    end
    % [DEC,AZM,RAD] = soradna(lat(1),lon(1),jday,yy);
    trise=rs(:,:,1);tset=rs(:,:,2); clear rs;
    daylen=trise; % initial value
    daylen(tset>trise)=tset(tset>trise)-trise(tset>trise);
    daylen(tset<trise)=24-trise(tset<trise)+tset(tset<trise); 
elseif sunrs_switch==0 % for MODIS ETM rs is derived  from following equations.
    [rs] = suncycle(lat(1),lon(1),dt);   % output rs in GMT time
    trise=rs(1);tset=rs(2); 
    if tset>trise
        daylen=tset-trise;
    elseif tset<trise  
        daylen=24-trise+tset;
    end
end
%##########################################################################
% compute surface solar radiation each hour
%##########################################################################
[R0,hsun]=solarR0(lon0,lon,lat,yy,mm,dd,hh,mn,ss); % R0: horizontal extraterrestial solar insolation,hsun:the height of the sun (rad) 
%R0=R0./max(1.0e-6,sin(hsun));  % none-horizontal extraterrestial

% [trise,tset,daylen]=shinetime(lon0,lon,lat,yy,mm,dd);   %CALL shinetime(trise,tset,daylen,lon0,lon,lat,yy,mm,dd)
if max(trise(:))<max(tset(:))&&(hh<=max(trise(:)))||(hh>=max(tset(:))) % for nighttime
    Rsfc = zeros(size(DEM));
    Rb= zeros(size(DEM));
    Rd= zeros(size(DEM));
    Rr= zeros(size(DEM));
    taub=zeros(size(DEM));
    taud=zeros(size(DEM));
    taur=zeros(size(DEM));
    return;
elseif max(trise(:))>max(tset(:))&&hh>=max(tset(:))&&hh<=max(trise(:))
    Rsfc = zeros(size(DEM));
    Rb= zeros(size(DEM));
    Rd= zeros(size(DEM));
    Rr= zeros(size(DEM));
    taub=zeros(size(DEM));
    taud=zeros(size(DEM));
    taur=zeros(size(DEM));    
    return;
end

if length(DEM)==1  % evaluate the code at point
    [Rsfc,Rb,Rd,Rr,taub,taud,taur]=Rad_hour(lon0,lon,lat,yy,mm,dd,hh,pa,ta,rh,albedo,shine);
    return
end

[taub,taud] = trans(yy,jday,lat,lon,hsun,pa,ta,rh); % taub: Beam radiation transmittance Eq. (9a); taud: Diffuse radiation transmittance	Eq. (9b)
taur        = 0.271+0.706*taub; % reflectance transmitivity,Kumar 1997, 

%#####################################################################
%     Compute slope and aspect,
%#####################################################################
if length(DEM)==1                           %single point
    slop=zeros(size(DEM));
    asp=zeros(size(DEM));
else
     [h,asp,slop] = hillshade(DEM,lon(1,:)*111*1000,lat(:,1)*111*1000); % lat*111*1000,transfer the latitude to meter distance,DEM is in meter unit
end
L      = lat;
L      = deg2rad(L);       % degree to radians conversion factor
sinL = sin(L);
cosL = cos(L);
tanL = tan(L);
sinSlop  = sin(slop);
cosSlop  = cos(slop);
cosSlop2 = cosSlop.*cosSlop;
sinSlop2 = sinSlop.*sinSlop;
sinAsp   = sin(asp);
cosAsp   = cos(asp);
term1 = ( sinL.*cosSlop - cosL.*sinSlop.*cosAsp);
term2 = ( cosL.*cosSlop + sinL.*sinSlop.*cosAsp);
term3 = sinSlop.*sinAsp;
%dS     = 23.45 * sin(2*pi* ( (284+jday)/tau_a ) );         % Solar declination, in degree, Kumar 1997, equation 4
dS     = 0.409 * sin(2*pi*((jday-93.5)/tau_a));             % Solar declination, in radians  Lond Di et al. 2010,eq. A2
% The solar declination varies from -23.44° at the (northern hemisphere) winter solstice, 
%  through 0° at the vernal equinox, to +23.44° at the summer solstice.
hsr    = real(acos(-tanL*tan(dS)));                        % angle at sunrise or sunset, in radians,  Kumar 1997, equation 6
                 

%########################################
% 消除地形的影响
%########################################
hh1      =  hh+mn/60+ss/3600;
hs       =trise;  %inital value
hs(hh1>trise)=hsr(hh1>trise)-(pi*(hh1-trise(hh1>trise))./daylen(hh1>trise));       % hour angle
hs(hh1<trise)=hsr(hh1<trise)-(pi*(hh1+(24-trise(hh1<trise)))./daylen(hh1<trise));  % hour angle

cos_i    = (sin(dS).*term1) + (cos(dS).*cos(hs).*term2) + (cos(dS).*term3.*sin(hs));
cos_i(cos_i<0)=0;    

cos_i(cos_i>1|cos_i<0)=NaN;  % solar zenith angle [0, pi/2]
Rb       = R0.*taub.*cos_i;                           % Kumar 1997, equation 15, 16
% Rb(Rb<0) = 0;
sinAlpha = sinL.*sin(dS)+cosL.*cos(dS).*cos(hs);      % Kumar 1997, equation 2,
sinAlpha(sinAlpha<=0|sinAlpha>1)=NaN;                 % altitude angle is negative when the sun drops below the horizon. pi>=alpha>=0
Rd       = R0.*taud   .* ((cos(slop)).^2)./ (2* sinAlpha);  % slope=0时，会使得 ((cos(slop)).^2)./ (2* sinAlpha) 的最非常大
Rr       = R0.*albedo.*taur .* ((sin(slop)).^2)./ (2* sinAlpha) ;
Rsfc     = Rb+Rd+Rr;

end

