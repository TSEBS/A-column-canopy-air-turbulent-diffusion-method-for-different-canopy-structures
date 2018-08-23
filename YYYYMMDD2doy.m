function [doy] = YYYYMMDD2doy(yy, mm, dd)
% YYYYMMDD2doy   Convert the year, month and day to day of year
% ########################################################################% 
% Zhaoze Deng, Zhaoze.Deng@gmail.com
% CPAC, Department of Atmospheric Sciences, Peking University
% 2008-07-01

noday=[31;28;31;30;31;30;31;31;30;31;30;31];
for i=1:12
    nday(1:length(yy),i)=noday(i);
end

nday(:,2) = eomday(yy, 2);

if dd>nday(mm)
    error('The input day is larger than the maximum days of the month.');
end

doy = dd;
for i=1:length(yy)
    mon = 1;
    while mon<mm(i)
        doy(i) = doy(i)+nday(i,mon);
        mon = mon+1;
    end
end


clear noday nday;