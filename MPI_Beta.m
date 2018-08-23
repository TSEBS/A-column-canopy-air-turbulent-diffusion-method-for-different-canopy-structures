function [beta]=MPI_Beta(jday,lat1,lon1,rh1)
% purpose: calculate the turbidity coef. according to the dataset of hess et al.complex(1998, bams);

%  input:
%  jday: julian day, 
%  lat1: local latitude, positive for northern hemisphere (deg.), length =1
%  lon1: longitude of each site (deg.), length =1, > 0 in East Hemesphere, < 0 in West Hemesphere
%  rh1: relative humidity (%),length =1

%  output:;
%  Angstrom turbidity coefficient

X=(-180:5:175);
Y=(90:-5:-90);
Z=[0,50,70,	80,	90,	95,	98,	99];

load betaw1;load betas1;
[y,x,z]=meshgrid(Y,X,Z);
if(jday>=60&&jday<=240)
    beta=griddata3(x,y,z,betas1,lon1,lat1,rh1); 
else
    beta=griddata3(x,y,z,betaw1,lon1,lat1,rh1);     
end

clear betaw1;
clear betas1;
return;


