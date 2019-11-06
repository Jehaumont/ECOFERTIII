function soilStateCilinderParams = store_cm(managementSettings,simulationSettings,soilStateCilinderParams)

%% FUNCTION INPUT
% management settings
nsol = managementSettings.nsol;

% simulation settings
print_node = simulationSettings.print_node;
print_time = simulationSettings.print_time;
simsol = simulationSettings.simSolFlag;
t = simulationSettings.t;

% soil state cilinder parameters
cm = soilStateCilinderParams.cm;
cm_profile = soilStateCilinderParams.cm_profile;
cm_node = soilStateCilinderParams.cm_node;
no_conv = soilStateCilinderParams.no_conv;


%% FUNCTION MAIN BODY
%Store concentration and temperature for output
if no_conv~=1	
	%wc_profile information;
    for j=1:nsol
        if simsol
            cm_node{1,j}(:,end+1)=(cm(print_node,j));
        end
    end
end

%% FUNCTION OUTPUT

% soil state cilinder parameters
soilStateCilinderParams.cm_node = cm_node;
soilStateCilinderParams.cm_profile = cm_profile;