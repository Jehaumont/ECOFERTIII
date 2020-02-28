function [climateState,managementSettings,simulationSettings,soilConsParams,soilStateCilinderParams]=...
        runwave(climateState,cropStateParams,cropConsParams,managementSettings,simulationSettings,soilConsParams,soilStateCilinderParams)

%% DOCUMENTATION

%% FUNCTION INPUT

sim_nitro = simulationSettings.simNitroFlag;
simplant = simulationSettings.simPlantFlag;
simsol = simulationSettings.simSolFlag;
simtemp = simulationSettings.simTempFlag;

%% FUNCTION MAIN BODY
% Evaluate bc at time t and t+dt, chooses right bc, checks if bc 
    % changed during previous iteration    
    soilStateCilinderParams = soil_boundary_conditions_new(climateState,simulationSettings,soilConsParams,soilStateCilinderParams);

    % Calculate ph
	soilStateCilinderParams = solve_flow(cropStateParams,simulationSettings,...
    managementSettings,soilConsParams,soilStateCilinderParams);
      
    %Calculate temperature
    if simtemp
        soilStateCilinderParams = solve_temperature(climateState,simulationSettings,soilConsParams,...
            soilStateCilinderParams);
    end
          
    %Calculate solute
    if simsol
        soilStateCilinderParams = wat_sol(managementSettings,simulationSettings,soilConsParams,soilStateCilinderParams);
        
        soilStateCilinderParams  = solve_solute(climateState,cropConsParams,cropStateParams,...
            managementSettings,simulationSettings,soilConsParams,soilStateCilinderParams);
        
       [soilConsParams,soilStateCilinderParams] = sol_intgr( managementSettings,...
           simulationSettings,soilConsParams,soilStateCilinderParams);
        if sim_nitro
             soilStateCilinderParams = nit_intg(simulationSettings,soilStateCilinderParams);
        end
        
    end    
    
soilStateCilinderParams =  wat_intgr(simulationSettings,soilConsParams,soilStateCilinderParams);

% Calculation of delta t
simulationSettings.t = simulationSettings.t + soilStateCilinderParams.dt;

soilStateCilinderParams = calc_dt_new(simulationSettings, soilStateCilinderParams); 
if simsol
    soilStateCilinderParams = calc_dt_sol(climateState,managementSettings,simulationSettings,soilConsParams,soilStateCilinderParams);
else
    soilStateCilinderParams.dt_sol = soilStateCilinderParams.dt;

end
soilStateCilinderParams.dt = min(soilStateCilinderParams.dt, soilStateCilinderParams.dt_sol);

% Evaluate bc
soilStateCilinderParams.bc(end+1,:) = ...
    [simulationSettings.t,...
    soilStateCilinderParams.boco_top_type,...
    soilStateCilinderParams.boco_top,...
    soilStateCilinderParams.boco_bot_type,...
    soilStateCilinderParams.boco_bot,...
    soilStateCilinderParams.pond,...
    soilStateCilinderParams.runoff,...
    soilStateCilinderParams.flxa1,...
    soilStateCilinderParams.dt,...
    soilStateCilinderParams.iter];

try
    soilStateCilinderParams.bc = soilStateCilinderParams.bc(end-1:end,:);
catch
    return
end

% Break when no convergence is reached
if soilStateCilinderParams.no_conv
    disp('No convergence');
	soilStateCilinderParams.noConvCounter = soilStateCilinderParams.noConvCounter+1;
    if soilStateCilinderParams.noConvCounter < 3
        return
    else
        soilStateCilinderParams.STOP = 1;
    end
end


