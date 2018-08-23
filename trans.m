function [taub,taud]=trans(yy,jday,lat,lon,hsun,pa,ta,rh)
%#######################################################################;
%     purpose:;
%     calculate the transmittance in clear skies;
%     Reference:Yang et 2006, Chen et al. 2012
%#######################################################################;
%     author: Xuelong Chen;

p0     = 101300. ;
%     calculate air mass;
mass   = 1./(sin(hsun)+0.15.*(57.3.*hsun+3.885).^(-1.253));
%   calculate turbidity;
[beta] = MPI_Beta(jday,lat(1),lon(1),rh);
%     calculate ozone thickness;
[loz]  = NASA_ozone(yy,jday,lat(1),lon(1));
%     calculate precipitable water amount;
water = 0.00493.*rh./ta.*exp(26.23-5416./ta);
%     calculate transmittance;
mb  = mass .* beta;
moz = mass .* loz;
mw  = mass .* water;
mc   = mass.*pa./p0;
taug = exp(-0.0117*mc.^0.3139);
taur = exp(-0.008735.*mc.*(0.547+0.014*mc-0.00038*mc.^2+0.0000046*mc.^3).^(-4.08));
tauw = min(1,0.909-0.036*log(mw));
tauoz= exp(-0.0365*moz.^0.7136);
taua = exp(-mb.*(0.677+0.1464*mb-0.00626*mb.^2).^(-1.3));

taub = max(0.0,taur.*taua.*tauoz.*taug.*tauw-0.013);            % Yang et 2006,eq.1a
taud = max(0,0.5 .*(tauoz.*taug.*tauw.*(1-taur.*taua)+0.013));  % Yang et 2006,eq.1b


return;