function simulationSettings = find_cropnumber (managementSettings, simulationSettings)
%% DOCUMENTATION
% This function find the number of the rotation which is being simulated at
% moment t. If no crop is being simulated the function returns NaN.

%% FUNCTION INPUT
%management settings
harvest_date = sort([managementSettings.hDateCauli;managementSettings.hDateLeek]); % set all harvest dates in chronological order
plant_date = sort([managementSettings.pDateCauli;managementSettings.pDateLeek]); % set all plant dates in chronological order

%simulation settings
simplant = simulationSettings.simPlantFlag; % boolean to turn off/on plant proces simulations
t = simulationSettings.t; % current time of simulation

%% FUNCTION MAIN BODY
isInGrowingSeason = t>=plant_date & t<= (harvest_date - 1); % check whether on the current simulation there should be a plant on the field

if any(isInGrowingSeason) % check if t is between planting and harvest date for any crop rotation
    ncrop = find(isInGrowingSeason); % find the number of the crop rotation for which t is inbetween harvest and planting
else
    ncrop = NaN;
end

if simplant && any(isInGrowingSeason)
    flag_double_sim = 1; 
else
    flag_double_sim = 0;    
end

%% FUNCTION OUTPUT
%simulation settings
simulationSettings.flag_double_sim = flag_double_sim; % boolean to indicate whether to simulate soil and plant OR only soil
simulationSettings.ncrop = ncrop; % number of rotation being simulated at moment t