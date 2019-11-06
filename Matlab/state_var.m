function [WC,kh,CH,rtex,EPRA]=state_var(cropStateParams,managementSettings,simulationSettings,...
                                        soilConsParams,soilStateCilinderParams,update_historic)
%% DOCUMENTATION
%Calculates state variables (moisture retention, hydraulic conductivity, differential 
%   moisture capacity and root extraction rates for a given pressure head)
%
%IN:
%   ph,soil_parameters,dt,update_historic,t,simplant,dx
%OUT:
% 	WC = volumetric soil moisture content
%   kh = soil unsaturated hydraulic conductivity (cm/min)
%   CH = differential moisture capacity (1/cm)
% 	rtex = root extraction rate (1/min)
%   EPRA=plant parameter
%CALL:
%	moist_ret,conduct,diff_moist_capac,RER
%CALLED BY:
%   solve_flow.m
%---------------------------
% M. Vanclooster 2/2/2000


%% FUNCTION INPUT
% management settings
harvest_date = sort([managementSettings.hDateCauli;managementSettings.hDateLeek]);
plant_date = sort([managementSettings.pDateCauli;managementSettings.pDateLeek]);

% simulation settings
t = simulationSettings.t;
dx = simulationSettings.dx;
simplant = simulationSettings.simPlantFlag;
units = simulationSettings.units;

% soil state cilinder params
epa = soilStateCilinderParams.epa;
ph = soilStateCilinderParams.ph;
soil_parameters = [soilConsParams.wcr,soilConsParams.wcs,soilConsParams.alfa,...
                  soilConsParams.N,soilConsParams.ks,soilConsParams.lambda,soilConsParams.alfa_r];

% crop state params
drz = cropStateParams.TotRootDepth(end);
fraction_plant = cropStateParams.fraction_plant(end);

%% FUNCTION MAIN BODY
WC = moist_ret(ph, soilStateCilinderParams, soilConsParams, update_historic);
kh=conduct(ph,soil_parameters);
CH=diff_moist_capac(ph,soil_parameters);
if simplant==1
   [rtex,EPRA] = RER(t,ph,simplant,dx,units,plant_date,harvest_date,drz,epa,fraction_plant);
else
   rtex=0.*ph';
   EPRA=0.*ph';
end
rtex=rtex';

