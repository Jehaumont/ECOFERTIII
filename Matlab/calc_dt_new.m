function soilStateCilinderParams =calc_dt_new(simulationSettings,soilStateCilinderParams)

%% DOCUMENTATION
% Calculates the optimal time step
%
%
% IN : 
%	dt : the actual time step (min)
%	iter : the number of "solve_flow" iterations
%   dt_min,dt_max
% OUT: 
% 	dt = the optimized time step (min)
% CALL:
%	none
% CALLED BY:
%   wavemat101.m
%-----------------------
% M. Vanclooster 2/2/2000

%% FUNCTION INPUT.
% simulation settings
dt_min = simulationSettings.dt_min;
dt_max = simulationSettings.dt_max;
t = simulationSettings.t;

% soil state cilinder parameters
BOUNDARY_CONDITIONS_MATRIX = soilStateCilinderParams.BOUNDARY_CONDITIONS_MATRIX;
dt = soilStateCilinderParams.dt;
iter = soilStateCilinderParams.iter;

%% FUNCTION MAIN BODY
if (iter<=5)
   dt=dt*2;
elseif ((iter>5)&&(iter<=10))
   dt=dt*1.6;
elseif ((iter>10)&&(iter<=20))
   dt=dt*1.3;
elseif (iter>25)
   dt=dt*0.7;
end

if dt<dt_min
   dt=dt_min;
elseif dt>dt_max
   dt=dt_max;
end

dt_start = dt;
% Adapt dt in terms of boundary conditions
i=max(find(BOUNDARY_CONDITIONS_MATRIX(:,1)<=t));        % at time t
i2=max(find(BOUNDARY_CONDITIONS_MATRIX(:,1)<=(t+dt)));  % at time t+dt
if i~=i2            % dt encompass two or more changes of BC
    if i~=i2-1      % dt encompass moe changes of BC
        i2=1+i;
	end
    dt=BOUNDARY_CONDITIONS_MATRIX(i2,1)-t; % reduction of the time step, this time step can be smaller than dtmin
    dt_changed_bc=1;
    if dt < dt_min  %if the new dt<dtmin, do the change in this iteration
        i=i2;
        dt_changed_bc=0;
        dt=dt_start;
    end
else
    dt_changed_bc=0;
end

%% FUNCTION OUTPUT
soilStateCilinderParams.dt = dt;
soilStateCilinderParams.dt_changed_bc = dt_changed_bc;
