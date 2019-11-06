function soilStateCilinderParams = store_solute_balance(managementSettings,simulationSettings,soilStateCilinderParams)

%% DOCUMENTATION
%Store the info calculated in sol_intgr in a variable named solute_balance.
%
% %IN:
% t: time
% tflsol : cumulative solute inflow (M L^-2)
% rleasa: cumulative netto flux at bottom of the profile (M L^-2)
% dsol: change in solute content for the entire profile (M L^-2)
% cberr_sol: error on the balance (M L^-2)
%tcsink : sink term (M L^-2)
% solute balance : summary of all above variables before
%
%OUT:
% solute balance : updated version of the solute balance

%% FUNCTION INPUT

% simulation settings
t = simulationSettings.t;

%  management settings
nsol = managementSettings.nsol;

% soil state cilinder parmeters
cberr_sol = soilStateCilinderParams.cberr_sol;
tflsol = soilStateCilinderParams.tflsol;
rleasa = soilStateCilinderParams.rleasa;
tcsink = soilStateCilinderParams.tcsink;
dsol = soilStateCilinderParams.dsol;
solute_balance = soilStateCilinderParams.solute_balance;
dt = soilStateCilinderParams.dt;

%% FUNCTION MAIN BODY
for i=1:nsol
    solute_balance{1,i}(end+1,:)= [t+dt tflsol(i) rleasa(i) dsol(i) cberr_sol(i) tcsink(i)];
end

%% FUNCTION OUTPUT
soilStateCilinderParams.solute_balance = solute_balance;

