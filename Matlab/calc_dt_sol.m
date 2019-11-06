function soilStateCilinderParams = calc_dt_sol(climateState,managementSettings,simulationSettings,soilConsParams,soilStateCilinderParams)
%% DOCUMENTATION
%Adapt dt in terms of solute boundary conditions
% IN:
% t : current time
% nsol : number of solutes
% dt: timestep
% applic_boolean: boolean indicating end of a solute application (1 = true, 0 = false)
% dt_min: minimum allowed timestep
% units: units of Mass, Time and Space
%OUT:
% dt_sol = timestep for the solute module
% applic_boolean: boolean indicating end of a solute application (1 = true, 0 = false)
%
%Joachim Vansteenkiste, 11/03/2010

%% FUNCTION INPUT
% climate state
climate = climateState.climateDaily;

% management settings
nsol = managementSettings.nsol;

% simulation settings
t = simulationSettings.t;
dt_min = simulationSettings.dt_min;
units = simulationSettings.units;

% soil constant parameters
ddepsol = soilConsParams.ddepsol;
wdepsol = soilConsParams.wdepsol;


% soil state cillinder parameters
applic_boolean = soilStateCilinderParams.applic_boolean;
conirsol = soilStateCilinderParams.conirsol;
csp = soilStateCilinderParams.csp;
dt = soilStateCilinderParams.dt;
dt_sol_count = soilStateCilinderParams.dt_sol_count;
fsol = soilStateCilinderParams.fsol;
reservoir = soilStateCilinderParams.reservoir;


%% FUNCTION MAIN BODY
%Find the current application
[conc,solute_applic,solboco_top_type,solboco_input_type]=In_solute_boundary(t,nsol);

if solboco_input_type == 2 %input is defined in M L-2
    
    %find current climatic conditions
    index = max(find(t>=climate(:,1)));
    current_rain = climate(index, 3);
    
    %Find current applications
    index = max(find(t>=fsol(:,1)));
    current_fert = fsol(index,:);
    index = max(find(t>=conirsol(:,1)));
    current_irr = conirsol(index,:);
    %Wet and dry deposition
    ressol=current_rain*dt.*wdepsol + ddepsol;
    %Inorganic fertiliser + Fertigation
    ressol = ressol + current_fert(2:end) + current_irr(2:end);
    cs = ressol + reservoir;
else
    cs = conc;
end


%If there is application decrease the timestep
if sum(cs) ~=0
    if strcmp(units{1,2},'day' )==1
        dt1 = 0.099;
    elseif strcmp(units{1,2},'hour' )==1
        dt1 = 0.0099*24;
    elseif strcmp(units{1,2},'min' )==1
        dt1 = 0.0099*24*60;
    elseif strcmp(units{1,2},'s' )==1
        dt1 = 0.0099*24*60*60;
    end
else
    dt1 = dt;
end

if solboco_input_type == 1
    %Adapt dt to the top BC
    i=max(find(solute_applic(:,1)<=t));        % at time t
    i2=max(find(solute_applic(:,1)<=(t+dt)));  % at time t+dt
    if i~=i2            % dt encompass two or more changes of BC
        if i~=i2-1      % dt encompass moe changes of BC
            i2=1+i;
        end
        dt2=solute_applic(i2,1)-t; % reduction of the time step, this time step can be smaller than dtmin
    else
        dt2 = dt;
    end
else
    dt2 =dt;
end

previous_time_step = sum(csp);
next_time_step = sum(cs);
if applic_boolean ==0 && isnan(previous_time_step)==0
    if previous_time_step ~=0 && next_time_step ==0 %application is finished, smaller time steps are required
        applic_boolean = 1;
        dt_sol_count = 0;
        dt3 = dt_min;
    else
        dt_sol_count = 0;
        dt3 = dt;
    end
else
    dt3 =dt;
end

if applic_boolean ==1
    dt_sol_count = dt_sol_count +1;
    dt3 = dt_min;
    if dt_sol_count == 10
        applic_boolean = 0;
    end
end
dt_sol = min(dt1, dt2);
dt_sol = min(dt_sol, dt3);

%% FUNCTION OUTPUT

% soil state cilinder parameters
soilStateCilinderParams.applic_boolean = applic_boolean;
soilStateCilinderParams.dt_sol = dt_sol;
soilStateCilinderParams.dt_sol_count = dt_sol_count;
