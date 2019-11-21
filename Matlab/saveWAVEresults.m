function soilStateCilinderParams = saveWAVEresults(climateState, fileSettings,...
    soilStateCilinderParams, simulationSettings, type)

%% FUNCTION INPUT
%simulation settings
t = simulationSettings.t;
dx = simulationSettings.dx;
%climate
climate = climateState.climateDaily;

%soil state cilinder params
bot_flux = soilStateCilinderParams.bot_flux;
bot_inf  = soilStateCilinderParams.bot_inf;
cberr_sol = soilStateCilinderParams.cberr_sol;
cm = soilStateCilinderParams.cm;
cum_bot_flxs = soilStateCilinderParams.cum_bot_flxs;
cum_evap = soilStateCilinderParams.cum_evap;
cum_infiltr = soilStateCilinderParams.cum_infiltr;
cum_pot_transp = soilStateCilinderParams.cum_pot_transp;
cum_potential_surface_flux  = soilStateCilinderParams.cum_potential_surface_flux;
cum_sink_wat = soilStateCilinderParams.cum_sink_wat1;
cum_sink_wat2 = soilStateCilinderParams.cum_sink_wat2;
cum_top_flxs = soilStateCilinderParams.cum_top_flxs;
cum_trans = soilStateCilinderParams.cum_trans;
water_balance = soilStateCilinderParams.water_balance;
dsol = soilStateCilinderParams.dsol;
epa = soilStateCilinderParams.epa;
esa = soilStateCilinderParams.esa;
evap = soilStateCilinderParams.evap;
flxar = soilStateCilinderParams.flxar;
ph = soilStateCilinderParams.ph;
PressHSoil = soilStateCilinderParams.PressHSoil;
potential_surface_flux = soilStateCilinderParams.potential_surface_flux;
potential_transp  = soilStateCilinderParams.potential_transp;
rleasa = soilStateCilinderParams.rleasa;
sink_prof = soilStateCilinderParams.sink_prof;
soil_om = soilStateCilinderParams.soil_om;
stock_initial = soilStateCilinderParams.stock_initial;
tcsolo = soilStateCilinderParams.tcsolo;
tcsink = soilStateCilinderParams.tcsink;
temp = soilStateCilinderParams.temp;
tflsol = soilStateCilinderParams.tflsol;
top_flux = soilStateCilinderParams.top_flux;
top_inf = soilStateCilinderParams.top_inf;
trans = soilStateCilinderParams.trans;
wat_flxs = soilStateCilinderParams.wat_flxs;
WatSinkSoil_log = soilStateCilinderParams.WatSinkSoil_log;
water_storage = soilStateCilinderParams.water_storage;
WC = soilStateCilinderParams.WC;
wcmo = soilStateCilinderParams.wcmo;
WCSoil_log = soilStateCilinderParams.WCSoil_log;
CarbMan_log = soilStateCilinderParams.CarbMan_log;
NitMan_log = soilStateCilinderParams.NitMan_log;
CarbLit_log = soilStateCilinderParams.CarbLit_log;
NitLit_log = soilStateCilinderParams.NitLit_log;
CarbHum_log = soilStateCilinderParams.CarbHum_log;
NitHum_log = soilStateCilinderParams.NitHum_log;
CO2vol_log = soilStateCilinderParams.CO2vol_log;
UreumSoil_log = soilStateCilinderParams.UreumSoil_log;
NH4Soil_log = soilStateCilinderParams.NH4Soil_log;
NO3Soil_log = soilStateCilinderParams.NO3Soil_log;
UreumConc_log = soilStateCilinderParams.UreumConc_log;
NH4Conc_log = soilStateCilinderParams.NH4Conc_log;
NO3Conc_log = soilStateCilinderParams.NO3Conc_log;
UreumBalance_log = soilStateCilinderParams.UreumBalance_log;
NH4Balance_log = soilStateCilinderParams.NH4Balance_log;
NO3Balance_log = soilStateCilinderParams.NO3Balance_log;
cum_nit_sink = soilStateCilinderParams.cum_nit_sink;
N_reaction_balance = soilStateCilinderParams.N_reaction_balance;
TSoil = soilStateCilinderParams.TSoil;
cberr_corg = soilStateCilinderParams.cberr_corg;
cberr_norg = soilStateCilinderParams.cberr_norg;
ptden = soilStateCilinderParams.ptden;
pthyd = soilStateCilinderParams.pthyd;
ptmin = soilStateCilinderParams.ptmin;
ptnit = soilStateCilinderParams.ptnit;
ptscorg = soilStateCilinderParams.ptscorg;
ptsnorg  = soilStateCilinderParams.ptsnorg; 
ptup = soilStateCilinderParams.ptup;
ptvol  = soilStateCilinderParams.ptvol;
tflcorg = soilStateCilinderParams.tflcorg;
tflnorg = soilStateCilinderParams.tflnorg;


%% FUNCTION MAIN BODY

% STORE INFORMATION IN TEMPORARY VARIABLES

% soil water content information
WCSoil_log(:,end+1) = WC.';
PressHSoil(:,end+1) = ph.';
WatSinkSoil_log(:,end+1) = cum_sink_wat2.';

%Cumulative bottom flux
bot_flux(end+1) = cum_bot_flxs;

%Instanate bottom flux
bot_inf(end+1) = wat_flxs(end);

%Cumulative top flux = cumulative actual surface flux
top_flux(end+1)= -cum_top_flxs;

%Instantanate top flux = Actual surface flux
top_inf(end+1)=  -wat_flxs(1);

%Potential surface flux
potential_surface_flux(end+1) = -flxar;

%Cumulative potential surface flux
cum_potential_surface_flux(end+1) = cum_potential_surface_flux(end) + flxar;

%Potential transpiration
potential_transp(end+1) = epa;

%Cumulative potential transpiration
cum_pot_transp(end+1) = cum_pot_transp(end) + epa;

%Actual transpiration;
trans(end+1) = sink_prof(end)*dx;

%cumulative actual transpiration
cum_trans(end+1) = sum(cum_sink_wat)*dx;

%Potential evaporation
evap(end+1) = esa;

i = find(climate(:,1)<=t, 1, 'last');
rinf = (climate(i,3) + climate(i,4));

if wat_flxs(1) > 0 && esa ~= 0
    
    cum_evap(end+1) = cum_evap(end) + abs(wat_flxs(1));
    cum_infiltr(end+1) = cum_infiltr(end);
    
    if rinf>0
        cum_evap(end) = cum_evap(end) + rinf;
        cum_infiltr(end) = cum_infiltr(end) + rinf;
    end
else
    
    cum_infiltr(end+1) = cum_infiltr(end) + abs(wat_flxs(1)) + esa;
    cum_evap(end+1) = cum_evap(end) + esa;
end

water_storage(end+1) = sum(WC.*dx);

water_balance(:,end+1)=[t top_flux(end), bot_flux(end), cum_trans(end), sum(wcmo)-stock_initial]';

% nitrogen

CarbMan_log(:,end+1) = soil_om(:,1);

NitMan_log(:,end+1) =  soil_om(:,2);

CarbLit_log(:,end+1) =  soil_om(:,3);

NitLit_log(:,end+1) =  soil_om(:,4);

CarbHum_log(:,end+1) =  soil_om(:,5);

NitHum_log(:,end+1) =  soil_om(:,6);

CO2vol_log(:,end+1) =   soil_om(:,7);

UreumSoil_log(:,end+1) = tcsolo(:,1);

UreumConc_log(:,end+1) = cm(:,1);

UreumBalance_log(:,end+1) = [t tflsol(1) rleasa(1) dsol(1) cberr_sol(1) tcsink(1)]';

NH4Soil_log(:,end+1) = tcsolo(:,2);

NH4Conc_log(:,end+1) = cm(:,2);

NH4Balance_log(:,end+1) = [t tflsol(2) rleasa(2) dsol(2) cberr_sol(2) tcsink(2)]';

NO3Soil_log(:,end+1) = tcsolo(:,3);

NO3Conc_log(:,end+1) = cm(:,3);

NO3Balance_log(:,end+1) = [t tflsol(3) rleasa(3) dsol(3) cberr_sol(3) tcsink(3)]';

cum_nit_sink(:,end+1) = [ptsnorg, ptscorg, ptup', ptmin', pthyd,ptnit,ptvol,ptden, tflnorg, tflcorg, cberr_norg, cberr_corg]';

N_reaction_balance(:,end+1) = [t ptup(1) ptup(2) ptup(3) ptmin(1) ptmin(2) ptmin(3) pthyd ptnit ptvol ptden]';

% temperature
TSoil(:,end+1) = temp;

%% SAVE VARIABLE TO STRUCTURE
soilStateCilinderParams.WCSoil_log = WCSoil_log;
soilStateCilinderParams.PressHSoil = PressHSoil;
soilStateCilinderParams.WatSinkSoil_log = WatSinkSoil_log;
soilStateCilinderParams.bot_flux = bot_flux;
soilStateCilinderParams.bot_inf = bot_inf;
soilStateCilinderParams.top_flux = top_flux;
soilStateCilinderParams.top_inf = top_inf;
soilStateCilinderParams.potential_surface_flux = potential_surface_flux;
soilStateCilinderParams.cum_potential_surface_flux = cum_potential_surface_flux;
soilStateCilinderParams.potential_transp = potential_transp;
soilStateCilinderParams.cum_pot_transp = cum_pot_transp;
soilStateCilinderParams.trans = trans;
soilStateCilinderParams.cum_trans = cum_trans;
soilStateCilinderParams.evap = evap;
soilStateCilinderParams.cum_evap = cum_evap;
soilStateCilinderParams.cum_infiltr = cum_infiltr;
soilStateCilinderParams.water_storage = water_storage;
soilStateCilinderParams.water_balance = water_balance;
soilStateCilinderParams.CarbMan_log = CarbMan_log;
soilStateCilinderParams.NitMan_log = NitMan_log;
soilStateCilinderParams.CarbLit_log = CarbLit_log;
soilStateCilinderParams.NitLit_log = NitLit_log;
soilStateCilinderParams.CarbHum_log = CarbHum_log;
soilStateCilinderParams.NitHum = NitHum_log;
soilStateCilinderParams.CO2vol_log = CO2vol_log;
soilStateCilinderParams.UreumSoil_log = UreumSoil_log;
soilStateCilinderParams.UreumConc_log = UreumConc_log;
soilStateCilinderParams.UreumBalance_log = UreumBalance_log;
soilStateCilinderParams.NH4Soil_log = NH4Soil_log;
soilStateCilinderParams.NH4Conc_log = NH4Conc_log; 
soilStateCilinderParams.NH4Balance_log = NH4Balance_log;
soilStateCilinderParams.NO3Soil_log = NO3Soil_log;
soilStateCilinderParams.NO3Conc_log = NO3Conc_log;
soilStateCilinderParams.NO3Balance_log = NO3Balance_log;
soilStateCilinderParams.cum_nit_sink = cum_nit_sink;
soilStateCilinderParams.N_reaction_balance = N_reaction_balance;
soilStateCilinderParams.TSoil = TSoil;

%% SAVE results to external storage
% set names of fields to save
fieldNames = {'JulianDay', 'WCSoil_log', 'PressHSoil', 'WatSinkSoil_log', 'bot_flux',...
              'bot_inf', 'top_flux', 'top_inf', 'potential_surface_flux',...
              'cum_potential_surface_flux', 'potential_transp',...
              'cum_pot_transp','trans', 'cum_trans', 'evap', 'cum_evap', 'cum_infiltr',...
              'water_storage', 'CarbMan_log', 'NitMan_log', 'CarbLit_log', 'NitLit_log',...
              'CarbHum_log', 'NitHum', 'CO2vol_log', 'UreumSoil_log', 'UreumConc_log',...
              'UreumBalance_log', 'NH4Soil_log', 'NH4Conc_log', 'NH4Balance_log',...
              'NO3Soil_log', 'NO3Conc_log', 'NO3Balance_log', 'cum_nit_sink',...
              'N_reaction_balance', 'TSoil','extra_sol', 'extra_water', ...
              'water_balance'};           

if simulationSettings.t == simulationSettings.tmax | soilStateCilinderParams.STOP
    soilStateCilinderParams = dataLog(soilStateCilinderParams, type,...
                                      fileSettings, 'force', fieldNames);
elseif simulationSettings.t < simulationSettings.tmax
    soilStateCilinderParams = dataLog(soilStateCilinderParams, type,...
                                     fileSettings, false, fieldNames);
end




