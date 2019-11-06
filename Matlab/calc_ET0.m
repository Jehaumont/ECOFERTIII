function climateState =calc_ET0(climateConsParams,climateState,cropData)
%% FUNCTION INPUT

climate = climateState.climate; % climate data [year,day,hour,radiation,temperature,relative humidity]

elev = climateConsParams.alt; % altitude
k = climateConsParams.karman; % von Karman constant
lat = climateConsParams.latitude; % latitude
long = climateConsParams.longitude; % longitude
sigma = climateConsParams.sigma; % stefan-boltzman constant
TZlong = climateConsParams.TZLong; % longitude of time zone
z_u = climateConsParams.z_u; % height of windspeed measurement
z_t = climateConsParams.z_t; % height of temperature and humidity measurement

% load crop related parameters
%LAI = cropData.LAIcov; % Leaf area index
%albedo = cropData.albedo; % albedo
%croph = cropData.cropHeight; % crop height in m
rlday = cropData.rlday;% bulk stomatal resistance during day for a well-illuminated leaf (s m-1)
rlnight = cropData.rlnight;% bulk stomatal resistance during night for a well-illuminated leaf (s m-1)

%Extract data out of climate data
time = jday2matlab(climate(:,2), climate(:,1))+ climate(:,3)/24;
Tair = climate(:,5);
RH = climate(:,6);
SlrMJm2 = climate(:,4)*60*60*10^-6;
u = RH./RH.*2;%Climatic file needs to be updated to also include wind speed.

% Specify crop variables (can be vectors for time series)
albedo = 0.23;
croph=0.12; % crop height (m); 0.12m is crop height assumed for grass in FAO56 
LAI = 24*croph;
% From Allen et al. (1998) p. 21:
% rl=100 s/m daily value for well-watered grass; this leads to rs=100/(0.5*24*0.12)=70 s/m as standard rs for grass in FAO56
% The bulk stomatal resistance, rl, is the average resistance of an individual leaf. This
% resistance is crop specific and differs among crop varieties and crop management. It usually
% increases as the crop ages and begins to ripen. There is, however, a lack of consolidated
% information on changes in rl over time for the different crops. The information available in
% the literature on stomatal conductance or resistance is often oriented toward physiological or
% ecophysiological studies.
% The stomatal resistance, rl, is influenced by climate and by water availability. However,
% influences vary from one crop to another and different varieties can be affected differently.
% The resistance increases when the crop is water stressed and the soil water availability limits
% crop evapotranspiration. Some studies indicate that stomatal resistance is influenced to some
% extent by radiation intensity, temperature, and vapour pressure deficit.

%% FUNCTION MAIN BODY

% Calculate extra-terrestrial solar radiation with matlab function Ra 
% (ASCE equations)
interval = time(2)-time(1); %interval between climate measurements in days
EndTime=time;
StartTime=time-interval; 
[Ra1,julday,beta,N]=Ra(lat,long, TZlong, StartTime, EndTime);

%Change Boltzmann constant to the current time interval
sigma = sigma*(interval*24); 

% Calculate clear-sky solar radiation
atmfactor=(0.75+2.E-5*elev);
Rso=atmfactor.*Ra1; % ASCE Eq. 47 p; 37

% Define function to calculate saturation vapour pressure (kPa)
calc_esat=@(T) 0.6108.*exp((17.27.*T)./(T+237.3)); % ASCE Eq. 37 p.29 (OK)

% calculate saturation and actual vapour pressure (kPa) from temperature and RH
esat=calc_esat(Tair);
eact=RH./100.*esat;

% Calculation of Rn 
% ******************

% Calculate ratio Rs:Rso
ntimes=size(time); % number of records in weather data
Rs2Rso=SlrMJm2./Rso; % Calculate ratio Rs:Rso

ifirst=find(gt(beta,0.3),1,'first'); % look for first period when sun angle is above 15 degrees
lastRs2Rso=SlrMJm2(ifirst)/Rso(ifirst);  % Initialize lastRs2Rso with 
% Rs:Rso ratio from that first period when sun angle is above 15 degrees

% Now check if angle of the sun above the horizon 
% at the midpoint of the time period is above 15° or 0.3 radians.
% If it is not, take the Rs:Rso ratio from the end of the previous day
for ii=1:ntimes
  if beta(ii)>0.3 % daytime
      lastRs2Rso=SlrMJm2(ii)/Rso(ii); % store value for night time 
  else             % night time
      Rs2Rso(ii)=lastRs2Rso; % take Rs:Rso ratio from the end of the previous day
  end
end

Rs2Rso=min(1.0,max(0.3,Rs2Rso)); % Limit the Rs:Rso value to the interval [0.3,1.0]
fcloud = 1.35.*Rs2Rso-0.35; % cloudiness function [dimensionless] (limited [0.05,1.0])
%fcloud = 1.35.*SlrMJm2./Rso-0.35;

% Calculate Rnl with Brunt equation (ASCE Eq. 44)
% Rnl = net outgoing long-wave radiation (MJ m-2) 
Rnl=sigma.*fcloud.*(0.34-0.14.*sqrt(eact)).*((Tair+273.16).^4);

% Calculate Rn
Rn = SlrMJm2*(1-albedo)-Rnl;% Net radiation (MJm-2)

% Calculation of G, gamma and DELTA
% *********************************

% daytime is defined as time when Rn>0
daytime=gt(Rn,0.0); 

% Soil heat flux (MJm-2) 
G=(0.1*daytime+0.5*(1-daytime)).*Rn;%ASCE Eq. 65 a-b (p. 44)

% Site specific calculations
% Calculate average pressure from elevation in metres (FAO56 Eq. 7)
AveragePressure=101.3*((293-0.0065*elev)/293).^5.26; % (kPa) OK

% Calculate psychometric constant from average pressure at site(FAO56 Eq.8)
gamma=0.665E-3*AveragePressure; % gamma=psychometric constant (kPa/°C) OK

% function to calculate the slope of the 
% Saturation Vapor Pressure-Temperature Curve (?) (kPa/°C)
calc_DELTA=@(T) 2503*exp((17.27*T)./(T+237.3))./(T+237.3).^2; % ASCE Eq. 36 p.28 (OK)
DELTA=calc_DELTA(Tair);

% Calculation of ET0
% ******************
d=0.67*croph;% zero plane displacement height (m)
% aerodynamic resistance (s m-1) ASCE Eq. B.2 p. B.5 (=appendix B) valid for neutral stability conditions
ra=log((z_u-d)/(0.123*croph)).*log((z_t-d)/(0.0123*croph))./(k^2*u);
% surface resistance (s m-1) ASCE Eq. B.3 on p; B.6
rs=(rlday.*daytime+rlnight.*(1-daytime))./(0.5*LAI);% Note that rl can be different for daytime and nighttime

% Calculate length of time period in seconds from period in days 
dtsec=24*60*60*(EndTime-StartTime);  

% Calculate ET0 with Penman-Monteith Eq. (ASCE Eq. 1 p. 4) and using box 6
% FAO 56, pag 26
ET0=(0.408*DELTA.*(Rn-G) + gamma*0.622/0.287./(Tair+273.16).*(esat-eact).*dtsec./ra)./(DELTA+gamma*(1+rs./ra)); % ET0 in mm over interval.

[timeday ET0_mm_per_day] = convert_to_daily_sum(ET0,time);

%% FUNCTION OUTPUT
climateState.ET0_cm_per_day = ET0_mm_per_day/10;
climateState.ET0_time = timeday;

%% HELPER FUNCTIONS
function [timeday, daily_var] = convert_to_daily_sum(var,time)

%Convert data to daily data by summation of all data of the same day
%Missing values are not noted as NaN, but are skipped in this conversion,
%for the reason that NaN can not be used as a boundary condition is wavemat
daynumber = floor (time);
minday = min(daynumber);
maxday = max(daynumber);
i=minday;
linenumber =1;
while i <= maxday
    index = find(daynumber==i);
    if isempty(index)==0
        if size(index,1)>=23
            s = cumsum(var(index));
            daily_var(linenumber) = s(end);
            timeday(linenumber)=i;
            linenumber = linenumber +1;
        end
    end
    i=i+1;
end


