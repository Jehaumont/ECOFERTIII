function soilStateCilinderParams = calc_stock(simulationSettings,soilConsParams,soilStateCilinderParams,update_historic)
%% DOCUMENTATION
% Calculates the stock of the whole sample
%
%IN:
%	ph = pressure head (L)
%OUT:
%	stock = water storage in the profile (L)
%CALL:
%   moist_ret.m
%CALLED BY
%   wavemat.m
%-------------------
% M. Javaux 15/05/2000

%% FUNCTION INPUT
%simualation settings
dx = simulationSettings.dx;

%soil inner state params
ph = soilStateCilinderParams.ph;


%% FUNCTION MAIN BODY   

WC = moist_ret(ph, soilStateCilinderParams, soilConsParams, update_historic);
stock = sum(WC.*dx);

%% FUCNTION OUTPUT
%soil inner state parameters
soilStateCilinderParams.stock_initial = stock;

