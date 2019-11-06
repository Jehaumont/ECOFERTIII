function soilStateCilinderParams = nit_intg(simulationSettings,soilStateCilinderParams)

%% DOCUMENTATION
%IN
% ptsnorg : cumulative organic nitrogen sink
%ptscorg :  cumulative organic carbon sink
%ptup : profile total solute uptake
%cupt : cumulative plant uptake in each compartiment
%ptmin :  profile total mineralisation
%cmin : cumulative mineralisation of nitrogen in each compartiment
%pthyd : profile total hydrolysis
%chyd : cumulative hydrolysis in each compartiment
%ptnit : profile total nitrification
%cnit : cumulative nitrification in each compartiment
%ptvol : profile total volatilisation
%cvol : cumulative volatilisation in each compartiment
%ptden : profile total denitrification
%cden : cumulative denitrification in each compartiment
%tnorgs : initial organic nitrogen profile total
%tcorgs : initial organic carbon profile total
%soil_om : actual organic carbon and notrogen in the profile
%tflnorg
%tflcorg
%OUT
%cberr_norg : mass balance error of the organic nitrogen
%cberr_corg : mass balance error of the organic carbon
%terms of the mass balance 

%% FUNCTION INPUT
% simulation settings
dx = simulationSettings.dx;

% soil state cilinder parameters
carblit = soilStateCilinderParams.carblit;
carbman = soilStateCilinderParams.carbman;
cden = soilStateCilinderParams.cden;
chyd = soilStateCilinderParams.chyd;
cmin = soilStateCilinderParams.cmin;
cnit = soilStateCilinderParams.cnit;
cupt = soilStateCilinderParams.cupt;
cvol = soilStateCilinderParams.cvol;
deccorg = soilStateCilinderParams.deccorg;
decnorg = soilStateCilinderParams.decnorg;
deniti = soilStateCilinderParams.deniti;
dt = soilStateCilinderParams.dt;
denitm = soilStateCilinderParams.denitm;
hydro_ureai = soilStateCilinderParams.hydro_ureai;
hydro_uream = soilStateCilinderParams.hydro_uream;
mineri = soilStateCilinderParams.mineri;
minerm = soilStateCilinderParams.minerm;
nitrifi = soilStateCilinderParams.nitrifi;
nitrifm = soilStateCilinderParams.nitrifm;
om_balance = soilStateCilinderParams.om_balance;
ptmin =  soilStateCilinderParams.ptmin;
ptden =  soilStateCilinderParams.ptden;
pthyd =  soilStateCilinderParams.pthyd;
ptnit =  soilStateCilinderParams.ptnit;
ptscorg =  soilStateCilinderParams.ptscorg;
ptsnorg =  soilStateCilinderParams.ptsnorg;
ptvol =  soilStateCilinderParams.ptvol; 
ptup =  soilStateCilinderParams.ptup;
rnitlit =  soilStateCilinderParams.rnitlit;
rnitman = soilStateCilinderParams.rnitman;
soil_om =  soilStateCilinderParams.soil_om;
tflcorg = soilStateCilinderParams.tflcorg;
tflnorg = soilStateCilinderParams.tflnorg;
tcorgs = soilStateCilinderParams.tcorgs;
tnorgs = soilStateCilinderParams.tnorgs;
uptakei = soilStateCilinderParams.uptakei;
uptakem = soilStateCilinderParams.uptakem;
volati = soilStateCilinderParams.volati;
volatm = soilStateCilinderParams.volatm;

%% FUNCTION MAIN BODY

c_manp = soil_om(:,1);
n_manp = soil_om(:,2);
c_litp = soil_om(:,3);
n_litp = soil_om(:,4);
c_hump = soil_om(:,5);
n_hump = soil_om(:,6);
co2 = soil_om(:,7);
%initialisation
volume=dx;
%profile total
ptclit = sum(c_litp);  %%actual profile total
ptnlit = sum(n_litp);  %%actual profile total
ptcman = sum(c_manp);  %%actual profile total
ptnman = sum(n_manp);  %%actual profile total
ptchum = sum(c_hump);  %%actual profile total
ptnhum = sum(n_hump);  %%actual profile total
ptco2 = sum(co2);
ptsnorg = ptsnorg + sum(decnorg)*volume*dt;  %% cumulative profile total sink (ptsnorg conservé)
ptscorg = ptscorg + sum(deccorg)*volume*dt;  %% cumulative profile total sink (ptscorg conservé)

ptnitorg=ptnlit+ptnman+ptnhum;      %% stock actuel d'azote
ptcarborg=ptclit+ptcman+ptchum;      %%stock actial de carbone   

%profile sink

%ammoniuum uptake
ptup(2) = ptup(2)+ (sum(uptakem(:,2))+ sum(uptakei(:,2)))*dx*dt;  % profile total uptake
cupt(:,2) = cupt(:,2)+ (uptakem(:,2) + uptakei(:,2))*dx*dt;       % cumulative uptake
%ammonium mineralisation
ptmin(2) = ptmin(2) + (sum(minerm(:,2))+ sum(mineri(:,2)))*dx*dt; % profile total mineralisation
cmin(:,2) = cmin(:,2) + (minerm(:,2)+mineri(:,2))*dx*dt;           % cumulative mineralisation
%nitrate uptake
ptup(3)= ptup(3)+ (sum(uptakem(:,3))+ sum(uptakei(:,3)))*dx*dt;  % profile total uptake
cupt(:,3) = cupt(:,3)+ (uptakem(:,3) + uptakei(:,3))*dx*dt;       % cumulative uptake
% nitrate mineralisation
ptmin(3) = ptmin(3) + (sum(minerm(:,3))+ sum(mineri(:,3)))*dx*dt; % profile total mineralisation
cmin(:,3) = cmin(:,3) + (minerm(:,3)+mineri(:,3))*dx*dt;           % cumulative mineralisation

% hydrolise
pthyd = pthyd + (sum(hydro_uream) + sum(hydro_ureai))*dx*dt;  % profilr total
chyd = chyd + (hydro_uream + hydro_ureai)*dx*dt;

% nitrification
ptnit = ptnit + (sum(nitrifm) + sum(nitrifi))*dx*dt;
cnit = cnit + (nitrifm + nitrifi)*dx*dt;
% volatilisation
ptvol = ptvol + (sum(volatm) + sum(volatm))*dx*dt;
cvol = cvol + (volatm + volati)*dx*dt;
% denitrification
ptden = ptden + (sum(denitm) + sum(deniti))*dx*dt;
cden = cden + (denitm + deniti)*dx*dt;

% organic balance term
dnorg=ptnitorg - tnorgs; %%(variation stock = actual-initial)
dcorg=ptcarborg - tcorgs; %% (variation stock = actual-initial)

% organic mass balance error
tflnorg = tflnorg + rnitlit + rnitman; %Input (defined in in_om_applic)
tflcorg = tflcorg + carblit + carbman; %Input (defined in in_om_applic)

cberr_norg= dnorg - (tflnorg - ptsnorg); %% (error = variation stock -(entrées-sorties))
cberr_corg= dcorg - (tflcorg - ptscorg);

%% FUNCTION OUTPUT

soilStateCilinderParams.cberr_corg	=	cberr_corg;
soilStateCilinderParams.cberr_norg	=	cberr_norg;
soilStateCilinderParams.cden	=	cden;
soilStateCilinderParams.chyd 	=	chyd;
soilStateCilinderParams.cmin	=	cmin;
soilStateCilinderParams.cnit 	=	cnit;
soilStateCilinderParams.cupt 	=	cupt; 
soilStateCilinderParams.cvol 	=	cvol; 
soilStateCilinderParams.dcorg	=	dcorg;
soilStateCilinderParams.dnorg	=	dnorg;
soilStateCilinderParams.om_balance	=	om_balance;
soilStateCilinderParams.ptden	=	ptden;
soilStateCilinderParams.pthyd	=	pthyd;
soilStateCilinderParams.ptmin	=	ptmin;
soilStateCilinderParams.ptnit	=	ptnit;
soilStateCilinderParams.ptscorg	=	ptscorg;
soilStateCilinderParams.ptsnorg 	=	ptsnorg; 
soilStateCilinderParams.ptup	=	ptup;
soilStateCilinderParams.ptvol	=	ptvol;
soilStateCilinderParams.tflcorg	=	tflcorg;
soilStateCilinderParams.tflnorg	=	tflnorg;

