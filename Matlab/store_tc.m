function soilStateCilinderParams = store_tc(managementSettings,simulationSettings,soilStateCilinderParams)

%% FUNCTION INPUT

% simulation settings
print_node = simulationSettings.print_node;
print_time = simulationSettings.print_time;
simsol = simulationSettings.simSolFlag;
t = simulationSettings.t;

% management settings
nsol = managementSettings.nsol;

% soil state parameters
no_conv = soilStateCilinderParams.no_conv;
tcsolo = soilStateCilinderParams.tcsolo;
tc_node = soilStateCilinderParams.tc_node;
tc_profile = soilStateCilinderParams.tc_profile;

%% FUNCTION MAIN BODY
%Store concentration for output (in M L-2)
if no_conv~=1	
        if simsol
            tc_node{1,j}(:,end+1)=(tcsolo(print_node,j));
        end
end

%% FUNCTION OUTPUT

soilStateCilinderParams.tc_node = tc_node;
soilStateCilinderParams.tc_profile = tc_profile;