function soilConsParams = PTFweynants(simulationSettings,soilCommonStateParams,soilConsParams)

%% DOCUMENTATION
% Pedotransfer function to estimate parameters of van Genuchten-Mualem
% model (VGM) with restriction m=1-1/n.  Implementation of Table 6 of
% Weynants et al., 2009
% Weynants, M., Vereecken, H., Javaux, M., 2009. Revisiting Vereecken
% Pedotransfer Functions: Introducing a Closed-Form Hydraulic Model. 
% Vadose Zone Journal 8, 86-95.
%
% Input:
% sand clay as percent
% carbon as gC/kg
% Bulk density as g cm-3
%
% output:
% wcr and wcs in cm³/cm³
% alpha in 1/cm
% Ks in cm/d
% n and lambda are dimensionless
%
% Regression equations in Table 6 of Weynants et al. (2009)

%% FUNCTION INPUT
% simulation setting
dx = simulationSettings.dx; % soillayer thickness

% soil common state parameters
soil_om = soilCommonStateParams.soil_om; % total carbon content in each soil layer

% soil constant parameters
bulkd = soilConsParams.bulkDensity; % bulk density in each soil layer
clay = soilConsParams.percClay; % percentage clay in each soil layer
sand = soilConsParams.percSand; % percentage sand in each soil layer

%% FUNCTION MAIN BODY
carbon = sum(soil_om(:,1:2:5),2)./dx; %[gC/cm^3]
carbon = carbon./(bulkd/1000); % [gC/kg], convert bulk density from g/cm^3 to kg/cm^3 

wcr =       zeros(size(bulkd));
wcs =       0.6355 + 0.0013*clay - 0.1631*bulkd;
alpha =     exp(-4.3003 - 0.0097*clay + 0.0138*sand - 0.0992*carbon);
n =         (1 + exp(-1.0846 - 0.0236*clay - 0.0085*sand + 0.0001*sand.^2)).*ones(size(bulkd));
lambda =    (-1.8642 - 0.1317*clay + 0.0067*sand).*ones(size(bulkd));

% Ks is calibrated because of convergence problems when using PTF and large
% paramter range
% Ks = exp(1.9582 + 0.0308*sand - 0.6142*bulkd -0.1566*carbon);
%% FUNCTION OUTPUT

% soil constant params
% soilConsParams.wcr = wcr;
soilConsParams.wcs = wcs;
soilConsParams.alfa = alpha;
%soilConsParams.N(1:30) = n(1:30);
soilConsParams.lambda = lambda;
%soilConsParams.ks = Ks;

