function [Rsfc,radb,radd,radr,taub,taud,taur]=Rad_hour(lon0,lon,lat,yy,mm,dd,hh,pa,ta,rh,albedo,shine)
%  INPUT:
%      real    lon0       ! the longitude where the standard time is defined. 
%                         ! If input time is UTC, lon0 = 0
%                         ! If input time is BST, lon0 = 120
%                         ! If input time is JST, lon0 = 135
%      real    lon        ! longitude of each site (deg.)
%                         ! > 0 in East Hemesphere
%                         ! < 0 in West Hemesphere
%      real    lat        ! latitdue of each site (deg.)
%                         ! > 0 in North Hemesphere
%                         ! < 0 in South Hemesphere
%      integer	yy,mm,dd,hh !calendar year:month:day:hour defined at lon0
%      real    pa         ! surface air pressure (Pa)
%      real    ta         ! air temperature (K)
%      real    rh         ! relative humidity
%      real    shine      ! Hourly sunshine duration (hour)

%  OUTPUT:
%     real    Rsfc        ! Hourly all-sky solar irradiance (W m^-2)

%  Temporary:
%      integer jday       ! julian day
%      integer	mn,ss     !standard time minute:second in lon0
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
if isnan(yy)||isnan(mm)||isnan(dd)||isnan(hh)||isnan(pa)||isnan(ta)||isnan(rh)||isnan(shine) % for NaN value
   Rsfc = NaN; R0=NaN; radd=NaN;radb=NaN;radr=NaN;
   taub=NaN;taud=NaN;
   return;
end

kk  = 3;    % number of integrating intervals
Rbn = 120;  % Critical normal beam to define sunshine
%      real time(kk), Rbni(kk),Rbi(kk),Rdi(kk)
%      real tbgn,tend
%      integer i, hh0
[jday] = YYYYMMDD2doy(yy, mm, dd);
%     calculate sunrise and sunset in hours
[trise,tset,daylen]=shinetime(lon0,lon,lat,yy,mm,dd);   %CALL shinetime(trise,tset,daylen,lon0,lon,lat,yy,mm,dd)
if(hh<trise)||(hh-1>tset) % for nighttime
   Rsfc = 0;R0=0;delta=0;radd=0;radb=0;radr=0;
   taub=0;taud=0;
   return;
end

tbgn  = max(trise,(hh-1));
tend  = min(tset, (hh));
dt    = (tend-tbgn)/(kk-1) ;
%     Calculate mean radiation over tbgn and tend by integrating 
radb = 0;
radd = 0;
radr = 0;
for i=1:kk
         time(i) = tbgn + dt * (i-1);
         hh0  = floor(time(i));%int(time(i));
         mn   = floor(time(i) * 60  - hh0 * 60);  %int(time(i) * 60  - hh0 * 60) ;
         ss   = floor(time(i) * 3600- hh0 * 3600 - mn * 60); %int(time(i) * 3600- hh0 * 3600 - mn * 60) ;
         [R0,hsun]=solarR0(lon0,lon,lat,yy,mm,dd,hh0,mn,ss); %R0: horizontal extraterrestial solar insolation,hsun:the height of the sun (rad)      
         %CALL solarR0(R0,hsun,lon0,lon,lat,yy,mm,dd,hh0,mn,ss)
         [taub,taud]=trans(yy,jday,lat,lon,hsun,pa,ta,rh); %taub: Beam radiation transmittance Eq. (9a); taud: Diffuse radiation transmittance	Eq. (9b)
         %CALL  trans(jday,lat,lon,alt,hsun,pa,ta,rh,taub,taud)
         Rbni(i)= R0*taub / max(1.0e-6,sin(hsun)) ;
         Rbi(i) = R0*taub;
         Rdi(i) = R0*taud;
         radb   = radb + Rbi(i)/kk;
         radd   = radd + Rdi(i)/kk;
         % Xuelong`s method
         taur   = 0.271+0.706*taub; % reflectance transmitivity,Kumar 1997, 
         radr   = radr + R0.*albedo.*taur;
end
rad    = radb + radd +radr ; % clear-sky surface solar irradiance
%     The following calculate maximum possible sunshine duration
%     Bright sunshine duration is defined by WMO as the time  
%     during which the direct solar radiation normal to the 
%     sun direction exceeds the level of 120 W/m2
shinemax = 0;
for i=2:kk
         if(Rbni(i-1)>Rbn)&&(Rbni(i)>Rbn)
            shinemax = shinemax + dt;
         elseif(Rbni(i-1)>Rbn)&&(Rbni(i)<Rbn)                              
            shinemax = shinemax   + dt* (Rbni(i-1)-Rbn)/(Rbni(i-1)-Rbni(i));
         elseif(Rbni(i-1)<Rbn)&&(Rbni(i)>Rbn)                              
            shinemax = shinemax   + dt* (Rbni(i)-Rbn)/(Rbni(i)-Rbni(i-1));
         end
end
%     Calculate all-sky solar irradiance,
shine = min(shine, shinemax);
ratio = shine/max(0.001,shinemax);
if(ratio>0) 
         Rsfc = (0.4435+0.3976*ratio+0.1589*(ratio^2))* rad;  %calibrate cloud-related transmittance
else  
         Rsfc = 0.2560*rad;
end



