function [month day] = doy2MMDD(year, doy)
%doy2MMDD   Convert day of year to year, month and day
%   Detailed explanation goes here
% ########################################################################% 
% Zhaoze Deng, Zhaoze.Deng@gmail.com
% CPAC, Department of Atmospheric Sciences, Peking University
% 2008-07-01

nday=[31;28;31;30;31;30;31;31;30;31;30;31];
nday(2) = eomday(year, 2);

if doy>sum(nday)
    error('The input doy is larger than the total days of the year!!');
end

month = 1;
day = doy;
while day>nday(month)
    day = day-nday(month);
    month = month+1;
end
