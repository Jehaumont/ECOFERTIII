function soilStateCilinderParams = store_T(simulationSettings,soilStateCilinderParams)

%% FUNCTION INPUT
%simulation settings
print_node = simulationSettings.print_node;
print_time = simulationSettings.print_time;
t = simulationSettings.t;

% soil state parameters
no_conv = soilStateCilinderParams.no_conv;
temp = soilStateCilinderParams.temp;
temp_node = soilStateCilinderParams.temp_node;
temp_profile = soilStateCilinderParams.temp_profile;

%% FUNCTION MAIN BODY
%Store concentration and temperature for output
if no_conv~=1
    temp_node(:,end+1)=(temp(print_node));
end

%% FUNCTIOIN OUTPUT

% soil state params
soilStateCilinderParams.temp_node = temp_node;
soilStateCilinderParams.temp_profile = temp_profile;