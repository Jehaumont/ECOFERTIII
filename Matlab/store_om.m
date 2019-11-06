function soilStateCilinderParams = store_om(simulationSettings,soilStateCilinderParams)

%% FUNCTION INPUT
% simulation settins
t = simulationSettings.t;
print_time = simulationSettings.print_time;
print_node = simulationSettings.print_node;
sim_nitro = simulationSettings.simNitroFlag;

% soil state params
no_conv = soilStateCilinderParams.no_conv;
om_node = soilStateCilinderParams.om_node;
om_profile = soilStateCilinderParams.om_profile;
soil_om = soilStateCilinderParams.soil_om;

%% FUNCTION MAIN BODY
%Store organic carbon and nitrogen content for output
if no_conv~=1
    for j=1:7 %there are 6 organic carbon and nitrogen pools & CO2
        if sim_nitro
            om_node{1,j}(:,end+1)= (soil_om(print_node,j));
        end
    end
end

%% FUNCTION OUTPUT

% soil state parameters
soilStateCilinderParams.om_node = om_node;
soilStateCilinderParams.om_profile = om_profile;