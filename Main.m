clc
clear all
constants1
global c1;global c2;global c3;global Ct;global hs1;global Cd;global fh;global A1;global A2;global A3;global xi_m;
global Hi_PBL;global beta0;global zs2hc;global dth;global phiuc; global hcp;global stp;global hs1;global As;
phiuc   =   2;
beta0   =   0.28;
zs2hc   =   3.45;
dth     =   0.7;
fh      =   0;
phiuc   =   3;
Cd      =   0.2;        % Foliage drag coefficient
Ct      =   0.01;       % Heat transfer coefficient
c1      =   0.320;                                           % model constants (Massman 1997)
c2      =   0.264;                                           % model constants (Massman 1997)
c3      =   15.1;                                            % model constants (Massman 1997)
Hi_PBL  =   1000;
hs1     =   0.004; % momentum roughness parameter (0.009 ~ 0.024)(Su et al., 1997, IJRS, p.2105-2124.),default value is 0.0012
kbc     =   5;
kbs     =   5;
A2      =   -5;
As      =   0.5;
Sigma_SB =   5.678E-8;      % Stefan-Boltzmann's constant (W/m2/K4)
yy=2010;mm=4;dd=9;hh=2;mn=30;ss=0; % time, GMT
col= 760;lin=1394;          % line and column of station point in the satellite data map
Zref            =   20;     % unit m
Tref_K          =   278.35; % AWS temperature unit K
Uref            =   4.9;    % AWS wind speed unit m/s
Rhref           =   18;     % AWS relative humidity unit %
sunrs_switch    =   0;      % switch for sunrise and sunset,0 or 1
pathi='.\input\';
pathj=[pathi,'NDVI'];
load (pathj);
pathj=[pathi,'albedo'];
load (pathj);
pathj=[pathi,'T'];
load (pathj);
pathj=[pathi,'hc'];
load (pathj);
pathj=[pathi,'fc'];
load (pathj);
pathj=[pathi,'LAI'];
load (pathj);
pathj=[pathi,'emissivity'];
load (pathj);
pathj=[pathi,'DEM'];
load (pathj);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
asl      =    DEM(lin,col);
LST_K    =    T; clear T;                     % T in cent.
Tair     =    Tref_K+(DEM-asl)*(-0.0065);     % air temperature adjusted with altitude, lapse rate=-0.0065 K/m
Pref     =    101325 * exp(-DEM./8430);       % Equation 4 in Yang 2006.
[lin1,col1]         = size(NDVI);
Uref(1:lin1,1:col1) = Uref;
Rhref(1:lin1,1:col1)= Rhref;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dms  =  [28 0 42];   % south latitude, degree, minute, second
lats = dms2degrees(dms);
dms  = [86 9 5];    % left latitude, degree, minute, second
lons = dms2degrees(dms);

dms  = [29 0 54.9]; % north latitude, degree, minute, second
late = dms2degrees(dms);
dms  = [87 1 42.81];% right latitude, degree, minute, second
lone = dms2degrees(dms);

lon       = linspace(lons,lone,col1);
lat       = linspace(lats,late,lin1);
[lon,lat] = meshgrid(lon,lat);
lat       = flipud(lat); %latitude and longitude of each point
% derive the beam, diffuse and reflected radiation (Rb, Rd,Rr), and its transmittance (taub,taud,taur)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[SWd,Rb,Rd,Rr,taub,taud,taur]=DEMRad_hour1(lon,lat,yy,mm,dd,hh,mn,ss,Pref,Tair,Rhref(1),DEM,albedo,sunrs_switch,1); % size of Rh,
%%
Hi_PBL               =  1500;
Earef                = (Rhref./100).* VAPPREMATRIX(Tair-273.15) ;     % observed actual vapor pressure in Pa,
r                    = (0.622*Earef)./(Pref-Earef);                   % mixing ratio (Ref: http://amsglossary.allenpress.com/glossary/search?id=mixing-ratio1)
qaref                = r ./ (1+r);                                    % Specific humidity (kg/kg) (% Ref: http://amsglossary.allenpress.com/glossary/search?id=specific-humidity1)
e0                   = Earef/100;                                     % in hPa
LWd                  = 1.24*[e0./Tref_K].^(1/7) * Sigma_SB .* Tref_K.^4;        % Brutsaert, W.1975 On a derivable formula for long-wave radiation from clear skies.Water Resour. Res., 11,742-744
SWu                  = SWd.*albedo;
LWu                  = 5.67*10^(-8)*emissivity.*LST_K.^4;
Rn                   = SWd+LWd-SWu-LWu;
% Ground Heat Flux, Kustas et al 1993
C                    = single(0.34);    % amplitude of LAI, obtained through fit with SCOPE
beta                 = single(0.46); % extinction coefficient, obtained through fit with fit with SCOPE
G                    = Rn.*C.*exp( - beta* LAI); % Kustas, W.P., Daughtry, C.S.T. van Oevelen P.J., Analatytical Treatment of Relationships between Soil heat flux/net radiation and Vegetation Indices, Remote sensing of environment,46:319-330 (1993)
%%
fc=single(fc);NDVI=single(NDVI);LAI=single(LAI);hc=single(hc);Zref=single(Zref);Uref=single(Uref);Pref=single(Pref);Tair=single(Tair);LST_K=single(LST_K);qaref=single(qaref);Rn=single(Rn);Hi_PBL=single(Hi_PBL);Earef=single(Earef);G=single(G);
%% SEBS Algorithm
beta =  zeros(size(LAI)); % zeros represent leaf off;, one represent leaf on;
delh =  0.0001;
z01  =  0:delh:1;
nz   =  length(z01);
stp  =  nz;
z0m  =  LAI;
d0   =  LAI;
z0h  =  LAI;
[m,n]=  size(LAI); 
LCT(1:m,1:n)= 12; % land covers, Grass,
% water = 0; evergreen needleleaf forest = 1;evergreen broadleaf forest = 2;deciduous needleleaf forest = 3;
% deciduous broadleaf forest = 4;mixed forests = 5;closed shrubland = 6; open shrublands = 7;
% woody savannas = 8; savannas = 9; grasslands = 10; permanent wetlands = 11; croplands = 12;
% urban and built-up = 13; cropland/natural vegetation mosaic = 14; snow and ice = 15;
% barren or sparsely vegetated = 16; unclassified = 255;
for i=1:m
    for j=1:n
        [hcp]     =  leaflength(LCT(i,j));
        [A2, As, lmdars,Cd]     =  MassmanPa(LCT(i,j));
        [z0m(i,j), d0(i,j), z0h(i,j)]=  kb_1(LCT(i,j),fc(i,j),NDVI(i,j),LAI(i,j),hc(i,j),Zref,Uref(i,j),Pref(i,j),Tair(i,j),LST_K(i,j),qaref(i,j),kbc,kbs);
    end
end
tic
[Rn,G0,H,LE,ustar] =   EnergyBalance_mapRn(d0, z0m, z0h,LCT,hc, fc,LAI, ...,
    albedo,emissivity,SWd,LWd, LST_K,Hi_PBL, Zref, Tair, Uref, Earef,qaref, Pref, Pref,0);
toc
%%
h1                                  =   figure('Position',[100 100 1124 900]);
h11                                 =   axes('Fontsize',12,'Nextplot','add');
set(gcf,'color','white')
h1=subplot(2,2,1);
hd=imagesc(real(H));
colorbar
title('(a) Sensible heat flux(W/m^2)','Fontsize',12)
set(gca,'ytick',[])
set(gca,'xtick',[])
grid on
caxis([-50 300])
hold on;
plot(col,lin,'+','Markersize',20,'color','w','LineWidth',2)

h2=subplot(2,2,2);
hd=imagesc(real(LE));
colorbar
title('(b) Latent heat flux(W/m^2)','Fontsize',12)
set(gca,'ytick',[])
set(gca,'xtick',[])
grid on
hold on;
plot(col,lin,'+','Markersize',20,'color','w','LineWidth',2)
caxis([-20 150])

h3=subplot(2,2,3);
hd=imagesc(real(Rn));
colorbar
title('(c) Net radiation(W/m^2)','Fontsize',12)
set(gca,'ytick',[])
set(gca,'xtick',[])
grid on
hold on;
plot(col,lin,'+','Markersize',20,'color','w','LineWidth',2)
caxis([-100 600])

h4=subplot(2,2,4);
hd=imagesc(real(G0));
colorbar
title('(d) Ground heat flux(W/m^2)','Fontsize',12)
set(gca,'ytick',[])
set(gca,'xtick',[])
grid on
hold on;
plot(col,lin,'+','Markersize',20,'color','w','LineWidth',2)
caxis([-50 200])
%%
saveas(gcf,'.\output\fluxes.emf', 'emf');
%% daily or monthly ET calculation
Lv           =   (2.501-0.00234*(Tair-273.15))*1000000;   % Eq. 2.20 in Nicolas, unit J kg_1
Rn_daily= 200; % W/m^2
ET_daily =8640000*evap_fr.*Rn_daily./(Lv*1000);
hd=imagesc(ET_daily);
colorbar
title('(a) Evapotranspiration (mm/d)','Fontsize',12)
set(gca,'ytick',[])
set(gca,'xtick',[])
grid on
hold on;
plot(col,lin,'+','Markersize',20,'color','w','LineWidth',2)
% caxis([-50 200])
%%
saveas(gcf,'.\output\Evapotranspiration_daily.emf', 'emf');