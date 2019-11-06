function soilStateCilinderParams =  wat_intgr(simulationSettings,soilConsParams,soilStateCilinderParams)

%% DOCUMENTATION
% Calculates the profile mass balance
%%%%%%%%%%%%%%%%%%%%
% M. Vanclooster 17/2/2000
%
% IN: 
%	ph = the soil pressure head
% 	wat_flxs = the soil moisture flux
%  cum_top_flxs = cumulative water flux at the top 
%  cum_bot_flxs = cumulative water flux at the bottom 
%  cum_sink_wat= cumulative sink of water 
% 	stock_initial = the initial soil moisture storage
%	dx,dt,ncomp
% OUT: 
% 	Mass balance error
% CALL:
% 	calc_stock,root_extraction_rates.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% FUNCTION INPUT
%simulation settings
ncomp = simulationSettings.nComp;

% soil state cilinder parameters
cum_bot_flxs =	soilStateCilinderParams.cum_bot_flxs;
cum_sink_wat =	soilStateCilinderParams.cum_sink_wat1;
cum_sink_wat2 =	soilStateCilinderParams.cum_sink_wat2;
cum_top_flxs =	soilStateCilinderParams.cum_top_flxs ;
dt =	soilStateCilinderParams.dt;
ph =	soilStateCilinderParams.ph;
rtex =	soilStateCilinderParams.rtex;
stock_initial =	soilStateCilinderParams.stock_initial;
wat_flxs =	soilStateCilinderParams.wat_flxs;

%% FUNCTION MAIN BODY
% Actual storage (cm)
tmp.ph = ph;
tmp.dt = dt;
tmp = calc_stock(simulationSettings,soilConsParams,tmp,0);
stock_actual = tmp.stock_initial;

% Cumulative flux at the top of the soil profile
% = Cumulative actual surface flux
cum_top_flxs = cum_top_flxs+wat_flxs(1)*dt;

% Cumulative flux at the bottom of the soil profile
cum_bot_flxs = cum_bot_flxs+wat_flxs(ncomp+1)*dt;

% Cumulative actual transpiration via actual transpiration (sink_prof)
sink_prof=sum(rtex);
cum_sink_wat = cum_sink_wat + sum(rtex)*dt;
cum_sink_wat2 = cum_sink_wat2 + dt*rtex;

% Change of storage (cm)
delta_stock = stock_actual-stock_initial;

% Mass balance error (cm)
mass_balance_error = delta_stock + cum_top_flxs - cum_bot_flxs - cum_sink_wat;

%% FUNCTION OUTPUT
soilStateCilinderParams.delta_stock	= delta_stock; 
soilStateCilinderParams.mass_balance_error = mass_balance_error;
soilStateCilinderParams.cum_top_flxs = cum_top_flxs;
soilStateCilinderParams.cum_bot_flxs = cum_bot_flxs; 
soilStateCilinderParams.cum_sink_wat = cum_sink_wat; 
soilStateCilinderParams.cum_sink_wat2 =	cum_sink_wat2;
soilStateCilinderParams.sink_prof =	sink_prof;
