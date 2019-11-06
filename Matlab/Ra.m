function [Ra,julday,beta,N]=Ra(lat,long, TZlong, StartTime, EndTime)
% Ra(lat,long, TZlong, StartTime, EndTime)
% Function to calculate extraterrestrial radiation Ra during the hour 
% (or shorter) period for a single location;  
% Ra = extraterrestrial radiation during period (MJ m-2 hour-1)
% Lat=latitude of site in decimal degrees (positive in northern hemissphere,  
% negative in southern) e.g. 40°30'N >> lat=40.5;   15°20'S >> lat=-15.3333 
% Long = longitude of the solar radiation measurement site [expressed as
% positive degrees WEST of Greenwich, England]
% TZLong = longitude of the center of the local time zone 
% [expressed as positive degrees WEST of Greenwich, England]. In the United
% States, TZlong = 75, 90, 105 and 120° for the Eastern, Central, Rocky Mountain
% and Pacific time zones, respectively, and TZlong = 0° for Greenwich, 
% 345° for Paris (France), and 255° for Bangkok (Thailand)
% StartTime, EndTime = (Matlab) Datenumber at start and end of period
% StartTime, EndTime can be vectors (calculation for time series) 

MidTime=(StartTime+EndTime)/2; % Calculate datenumber at midpoint of period
[year month day hr minutes seconds]=datevec(MidTime); % convert datenumber to year, month, day, ...
julday=day-32+floor(275*month/9)+2*floor(3./(month+1))+floor(month/100-mod(year,4)/4+0.975); % Eq. 52; julian day number (1 .. 365 or 366)

% Equation numbers refer to equations in ASCE (2005) (2002 draft)
Gsc=4.94; % solar constant = 4.92 MJ m-2 h-1
dr=1+0.033*cos(2*pi/365*julday); % Eq. 50; inverse relative distance factor (squared) for the earth-sun (unitless)
delta=0.409*sin(2*pi/365*julday-1.39); % Eq. 51; delta = solar declination (radians)
phi=(pi/180)*lat; % phi = latitude (radians); Eq. 49
omegas=acos(-tan(phi)*tan(delta));% sunset hour angle, ?s; Eq. 59

% seasonal correction for solar time (Eq. 57-58)
b=2*pi*(julday-81)/364;
Sc=0.1645*sin(2*b)-0.1255*cos(b)-0.025*sin(b);

% Calculate length of period (in hours)
dt=(EndTime-StartTime).*24;
% Time in hours at midpoint of period (value between 0.0 and 24.00)
tmid=mod(MidTime,1).*24; 

% solar time angle (radians) at the midpoint of the period (Eq. 55)
omega=pi/12*((tmid+0.06667*(TZlong-long)+Sc)-12); 

% omega1 = solar time angle at beginning of period (radians)
omega1=omega-pi/24*dt; % Eqs. 53
omega1=max(-omegas,min(omega1,omegas)); % ensure that omega1 is within [-omegas, omegas]
% omega2 = % solar time angle at end of period (radians)
omega2=omega+pi/24*dt; % Eqs. 54
omega2=max(-omegas,min(omega2,omegas)); % ensure that omega2 is within [-omegas, omegas]

% Calculation of Ra = extraterrestrial radiation during the hour (or shorter) period [MJ m-2 hour-1]
Ra=max(0,12./pi.*Gsc.*dr.*((omega2-omega1).*sin(phi).*sin(delta)+cos(phi).*cos(delta).*(sin(omega2)-sin(omega1)))); 
% (Eq. 48, and Ra can not become negative)

% angle of the sun above the horizon at the midpoint of the time period
beta=asin(sin(phi).*sin(delta)+cos(phi).*cos(delta).*cos(omega));

N = 24.*omegas./pi;

