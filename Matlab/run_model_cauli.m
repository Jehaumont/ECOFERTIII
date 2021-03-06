function [cropStateParams,soilInnerStateParams,soilOuterStateParams] =...
                                            run_model_cauli(cauliConsParams,...
                                            climateConsParams,climateState,...
                                            cropStateParams,fileSettings,...
                                            managementSettings,...
                                            simulationSettings,soilCommonStateParams,...
                                            soilInnerStateParams,soilOuterStateParams)
%% DOCUMENTATION

%% FUNCTION INPUT
% climate data
ET0_time = climateState.ET0_time;
ET0_cm_per_day = climateState.ET0_cm_per_day;
RefClim = climateState.crop_climate;

% climate constant parameters
AtmP = climateConsParams.atmp;
Cp = climateConsParams.Cp; %specific heat of air
CO2L = climateConsParams.CO2L; %CO2 level

% management settings
if simulationSettings.calibrateFlag ==1
    PLM2 = 4;
else
    PLM2 = managementSettings.PLM2Cauli; % plant density [plants/M?]
end

plant_date = sort([managementSettings.pDateCauli;managementSettings.pDateLeek]);

% soil common parameters 
NSupply = soilCommonStateParams.NSupply; % the amount of nitrogen [g] that is available for the plant to take up

% simulation settings
columndepth = simulationSettings.columnDepth;
ncrop = simulationSettings.ncrop; % number of rotation that is currently on the field
NFAST = simulationSettings.Nfast; % iterations in the fast loop per day
t = simulationSettings.t;% time of simulation
date = datevec(t);

TLayer = simulationSettings.dx; % thickness of a soil layer [cm]
DAP = t-plant_date(ncrop)+1;
DTFAST = 1.0/NFAST;

% cauliflower constant parameters
albedo = cauliConsParams.albedo;
CK = cauliConsParams.Ck; % exponential decay factor for vapour pressure deficit
CritN = cauliConsParams.CritN; % critical nitrogen concentration on cauliflower shoot
CurdFrac = cauliConsParams.CurdFrac; % fractional allocation of the assimilated dry matter to the curd related to the DVS [fraction, DVS]
CurdNc = cauliConsParams.CurdNc; % critical nitrogen concentration in the curd
DeadLeafFraction = cauliConsParams.DeadLeafFraction; % fraction of theleaves that die related to TeSum
RootLateralDistance = cauliConsParams.RootLateralDistance; %root lateral distance related to TeSum of RootTeSum
DVSeffect = cauliConsParams.DVSEffect; % reduction factor on maximal photosynthesis based on DVS [reduction factor,DVS]
Ilsolrad = cauliConsParams.Ilsolrad; % conversion factor of radiation to long wave radiation
Issolrad = cauliConsParams.Issolrad; % conversion factor of radiatino to short wave radiation
GREF = cauliConsParams.GREF; % growth efficiency g dry matter per g assimilate
LeafCoverageConversion = cauliConsParams.LeafCoverageConversion; % %soil coverd by leaves based on DVS [%coverage,DVS]
LeafFrac = cauliConsParams.LeafFrac; % fractional allocation of assimilated dry matter to the leaves related to DVS [%,DVS]
PARsolrad = cauliConsParams.PARsolrad; % %photosynthetically active radiation
PPFDpar = cauliConsParams.PPFDpar;
Q10 = cauliConsParams.Q10; % increase of maintenance respiration when temperature increases 10?C
Qe = cauliConsParams.Qe; % PAR utilisation coefficient [mgCO2/J]
RMRL = cauliConsParams.RMRL; % respiration coefficient for the leaves [gCH20/(gDM*Day)]
RMRF = cauliConsParams.RMRF; % respiration coefficient for the curd/stem [gCH20/(gDM*Day)]
RootFrac = cauliConsParams.RootFrac; % fraction of the total dry mass that is allocated to the roots
RootNc = cauliConsParams.RootNc; % critical nitrogen concentration in the root [%DM]
RootTeSum = cauliConsParams.RootTeSum; % effective temperature sum for root development
SLAEffect = cauliConsParams.SLAEffect; % sla converstion to leaf area based on DVS [converstion factor, DVS]
SRL = cauliConsParams.SRL; % specific root length [cm/g]
StemFrac = cauliConsParams.StemFrac; % fractional allocation of assimilated dry matter to the stem related to DVS [%, DVS]
TAU = cauliConsParams.TAU; % leaf conductance for C02 transfer [m/s]
Tb = cauliConsParams.Tb; % base temperature [?C]
TempEffect = cauliConsParams.TempEffect; % temperature reduction on maximal photosynthetic efficiency [reduction factor, temperature]
TeSumCurdIni = cauliConsParams.TeSumCurdIni; % effective temperature sum after which a curd is initiated [?C day]
TeSumCurdMat = cauliConsParams.TeSumCurdMat; % effective temperature sum after initiation whereafter the curd is matured [?C day]
Threshold = cauliConsParams.Threshold; % threshold for minimal root length
virRd = cauliConsParams.virRd; % relative length (to the soil layers) for root density distribution calculation
VPDL = cauliConsParams.VPDL; % vapour pressure deficit limit [kPa]
XKK = cauliConsParams.XKK; % light extinsion coefficient
XM = cauliConsParams.XM; % leaf transmission coefficient
x_y_param = cauliConsParams.x_y_param; % parameter for the lateral distribution of the roots
z1_param = cauliConsParams.z1_param;
z2_param = cauliConsParams.z2_param;

% cauliflower state parameters

CumCurdNuptake = cropStateParams.CumCurdNuptake(end); % cumulative curd nitrogen uptake [g/plant]
CumEvaporation = cropStateParams.CumEvaporation(end); % cumulative evaporation [g/plant]
CumLeafNuptake = cropStateParams.CumLeafNuptake(end); % cumulative leaf nitrogen uptake [g/plant]
CumLeafStemNuptake = cropStateParams.CumLeafStemNuptake(end); % cumulative leaf and stem nitrogen uptake [g/plant]
CumRootNuptake = cropStateParams.CumRootNuptake(end); % cumulative root nitrogen uptake [g/plant]
CumStemNuptake = cropStateParams.CumStemNuptake(end); % cumulative stem nitrogen uptake [g/plang]
CumShootNuptake = cropStateParams.CumShootNuptake(end); % cumulative shoot nitrogen uptake [g/plant]
CumTranspiration = cropStateParams.CumTranspiration(end); % cumulative transpiration
CurdNDemand = cropStateParams.CurdNDemand(end); % daily curd potential nitrogen demand [g/plant]
DMCurd = cropStateParams.DMCurd(end); % daily curd dry matter assimilation [g/plant]
DMLeaf = cropStateParams.DMLeaf(end); % daily leaf dry matter assimilation [g/plant]
DMStem = cropStateParams.DMStem(end); % daily stem dry matter assimilation [g/plant]
DMShoot = cropStateParams.DMShoot(end);
DVS = cropStateParams.DVS(end); % development stage
LeafNDemand = cropStateParams.LeafNDemand(end); % daily leaf nitrogen demand [g/plant]
RootNDemand = cropStateParams.RootNDemand(end); % daily root nitrogen demand [g/plant]
ShootNDemand = cropStateParams.ShootNDemand(end); % daily nitrogen demand of the whole plant [g/plant]
StemNDemand = cropStateParams.StemNDemand(end); % daily stem nitrogen demand [g/plant]
TeSum = cropStateParams.TeSum(end); % cumulative effective temperature
TotDM = cropStateParams.TotDM(end);
TotDMCurd = cropStateParams.TotDMCurd(end); % total curd dry matter [g/plant]
TotDMDeadLeaf = cropStateParams.TotDMDeadLeaf(end); % total dry matter of dead leaves [g/plant]
TotDMGreenLeaf = cropStateParams.TotDMGreenLeaf(end); % total dry matter of green leaves [g/plant]
TotDMLeaf = cropStateParams.TotDMLeaf(end); % total dry matter of leaves [g/plant]
TotDMLeafStem = cropStateParams.TotDMLeafStem(end);
TotDMStem = cropStateParams.TotDMStem(end); % total dry matter of the stem [g/plant]
TotDMRoot = cropStateParams.TotDMRoot(end); % total dry matter of the roots [g/plant]
TotNDemandCurd = cropStateParams.TotNDemandCurd(end); % total curd nitrogen demand [g/plant]
TotNDemandLeaf = cropStateParams.TotNDemandLeaf(end); % total leaf nitrogen demand
TotNDemandLeafStem = cropStateParams.TotNDemandLeafStem(end); % total nitrogen demand
TotNDemandRoot = cropStateParams.TotNDemandRoot(end); % total root nitrogen demand [g/plant]
TotNDemandStem = cropStateParams.TotNDemandStem(end); % total stem nitrogen demand [g/plant]
TotDMShoot = cropStateParams.TotDMShoot(end); % total dry matter of the whole plant [g/plant]

if DAP == 1
    % calculate initial LAI and leaf area
    TotLeafArea = DMLeaf * SLAEffect(1,1)/10000;
    NSupply = ShootNDemand;
    NFac = cropStateParams.NFac(1);
else
    TotLeafArea = cropStateParams.TotLeafArea(end);

end

%% FUNCTION MAIN BODY
% Soil nitrogen uptake and distribution over the organs
if ShootNDemand <= NSupply
    ShootNuptake = ShootNDemand;
    Novershoot = NSupply - ShootNDemand;
else
    ShootNuptake = NSupply;
    Novershoot = 0;
end

if ShootNDemand ~=0
    ratiostem = StemNDemand/ShootNDemand;
    ratioleaf = LeafNDemand/ShootNDemand;
    ratioroot = RootNDemand/ShootNDemand;
    ratiocurd = CurdNDemand/ShootNDemand;
else
    ratioleaf = 0;
    ratiostem = 0;
    ratioroot = 0;
    ratiocurd = 0;
end

LeafNuptake = ShootNuptake*ratioleaf;
StemNuptake = ShootNuptake*ratiostem;
RootNuptake = ShootNuptake*ratioroot;
CurdNuptake = ShootNuptake*ratiocurd;

% integration of total nitrogen uptake and calculation of Nfac
if DAP ~= 1
    
    % integration of total plant nitrogen uptake
    CumLeafNuptake = CumLeafNuptake + LeafNuptake;
    CumStemNuptake = CumStemNuptake + StemNuptake;
    CumLeafStemNuptake = CumStemNuptake + CumLeafNuptake;
    CumShootNuptake = CumShootNuptake + ShootNuptake;
    CumRootNuptake = CumRootNuptake + RootNuptake;
    CumCurdNuptake = CumCurdNuptake + CurdNuptake; 
    
    % calculation of nitrogen correction factor based weighted average of
    % all organs
    
    PercNLeafStem = max(0, CumLeafStemNuptake/TotDMLeafStem*100);
    PercNCurd = max(0, CumCurdNuptake/TotDMCurd*100);
    PercNRoot = max(0, CumRootNuptake/TotDMRoot*100);
   
    NFac_LeafStem = PercNLeafStem/interp1(CritN(:,2), CritN(:,1), DVS);
    NFac_Curd = PercNCurd/CurdNc;
    NFac_Root = PercNRoot/RootNc;
    
    DMratio_LeafStem = (TotDMGreenLeaf + TotDMStem)/TotDMShoot;
    DMratio_Curd = TotDMCurd/TotDMShoot;
    DMratio_Root = TotDMRoot/TotDMShoot;
    
    NFac = NFac_LeafStem * DMratio_LeafStem + NFac_Curd * DMratio_Curd + ...
           NFac_Root * DMratio_Root;   
    
    NFac = min(NFac,1);
    NFac = max(NFac, 0.05);
end

% Initialize total transpiration and evaporation and aily temperature sum
Tsum = 0;
GP = 0;
RMAINT = 0;

% FAST LOOP (HOURLY LOOP)
for JF = 1:NFAST
    % CALCULATION OF CLIMATE CONDITIONS
    TFAST = (JF-1) * 24/NFAST; % hour of the day
    d = find(RefClim(:,1)==t & RefClim(:,2)==(TFAST+1)); % indicator to get climate data of moment TFAST

    if length(d)>1
        disp('more than 1 found!!');
    end
    if isempty(d)
        disp(' invalid index of climate data');
    end

    solrad = RefClim(d,3);
    TMPA = RefClim(d,4);
    RH = RefClim(d,5);

    %VPD calculation in kPa
    Ps = 0.61078 * exp((17.2694 * TMPA)/(TMPA+237.3));

    %Calculate actual vapor pressure PV
    Pv = RH/100 * Ps;
    VPD = Ps-Pv;
    PAR = solrad * PARsolrad;
    PPFD = PAR * PPFDpar;
    
    % CALCULATION OF GROSS PHOTOSYNTHETIC PRODUCTION
    GPF = 0;
    
    % effect of CO2 on Pmax (from Gainesville)
    PMAX = TAU(1) * CO2L;
    if CO2L>1500
        PMAX = (TAU(1) * 1500) + (TAU(2) * (CO2L-1500));
    end
    
    % reduction of Pmax at extreme temperatures
    PMAX = PMAX * interp1(TempEffect(:,2),TempEffect(:,1),TMPA);
    
    %reduction of Pmax in relation to DVS, senescence
    PMAX = PMAX * interp1(DVSeffect(:,2),DVSeffect(:,1),DVS);
    
    %reduction of Pmax in relation to critical N concentration
    PMAX = PMAX * NFac;
    
    %reduction of XK with leaf area (from Olesen, 1997)
    if TotLeafArea < 0.235
        XK = 0.90 - 3.00 * TotLeafArea + 6.38 * TotLeafArea^2;
    else
        XK = XKK;
    end
    
    if PPFD >= 0.001
        % reduction of Pmax by VPD
        if VPD>=VPDL
            PMAX = PMAX * exp(CK*(VPD-VPDL));
        end
        % Acock's model
        TOP = (1-XM) * PMAX + Qe * XK * PPFD;
        BOT = (1-XM) * PMAX + Qe * XK * PPFD * exp(-XK * TotLeafArea * PLM2);
        GPF = (PMAX/XK) * log(TOP/BOT);
        
        % conversion from CO2 to CH2O (30/44 = 0.682)
        GPF = GPF * 0.682;
        % conversion of GPF from ?M/m2-s into g/m2-day   %?M/m2-s x 0.000044g/?M x 3600s/h x 24h/d = 3.8016 g/m2-day
        GPF = GPF * 3.8016/PLM2;
    end
    
    % CALCULATION OF MAINTENANCE RESPIRATION
    % Effect of temperature on maintenance respiration on hourly basis
    RMAINTF = (RMRL * (DMLeaf+DMStem)+RMRF * DMCurd) * (Q10^(0.1 * TMPA-2.0));
    
    % integration of variables on 24 hours
    GP = GP + GPF * DTFAST;
    RMAINT = RMAINT + RMAINTF * DTFAST;
    
    % calculate average day temperature
    Tsum = Tsum + TMPA;
end
% DAILY LOOP

% CALCULATE TOTAL NET INCREASE IN DRY MATTER
%#Transfert to Biomass, calculate Total DM in g per plant per dag
DM_daily = (GP - RMAINT) * GREF;

% CALCULATE DEVELOPMENTAL STAGE (DVS)
MeanDayTemp = Tsum/NFAST;
Td = MeanDayTemp - Tb;
%effective temperature set to zero if negative, non-decreasing Tsum
if Td > 0
    Te = Td;
else
    Te = 0;
end
%calculate Temperature sum
TeSum = TeSum + Te;

%DVS for leaves + stem
DVS_VegC = TeSum/TeSumCurdIni;
if DVS_VegC > 1
    DVS_Veg = 1;
    TeSumC = TeSum - TeSumCurdIni;
    DVS_Curd = TeSumC/TeSumCurdMat;
else
    DVS_Veg = DVS_VegC;
    DVS_Curd = 0;
end


%Overall DVS
DVS = DVS_Veg + DVS_Curd;
DVS = min(DVS,2);

% COMPARTIMENTATION OF DRY MATTER
%Growth rates g per plant per day
DMLeaf = DM_daily * interp1(LeafFrac(:,2),LeafFrac(:,1),DVS)/(1 + RootFrac);
DMStem = DM_daily * interp1(StemFrac(:,2),StemFrac(:,1),DVS)/(1 + RootFrac);
DMCurd = DM_daily * interp1(CurdFrac(:,2),CurdFrac(:,1),DVS)/(1 + RootFrac);

%Calculate DM of green and dead leaves
Kdl = interp1(DeadLeafFraction(:,2),DeadLeafFraction(:,1),Te);
if DVS > 1.2
    DMDeadLeaf = DMLeaf*(exp(Kdl)-1);
else
    DMDeadLeaf = 0;
end
if DVS > 1.2
    DMGreenLeaf = DMLeaf-DMDeadLeaf;
else
    DMGreenLeaf = DMLeaf;
end

%compartimentation to roots, 9 % of shoot and constant
DMShoot = DMGreenLeaf - DMDeadLeaf + DMStem + DMCurd;
DMRoot = DMShoot * RootFrac;

%calculate leaf area out of SLA-DVS curve
LeafArea = DMGreenLeaf * interp1(SLAEffect(:,2),SLAEffect(:,1),DVS)/10000; % [m^2]
% INTEGRATION OF RESULTS
TotDM = TotDM  + DM_daily;
TotDMLeaf = TotDMLeaf + DMLeaf;
TotDMStem = TotDMStem + DMStem;
TotDMCurd = TotDMCurd + DMCurd;
TotDMShoot = TotDMShoot + DMShoot;
TotDMRoot = TotDMRoot + DMRoot;
TotDMGreenLeaf = TotDMGreenLeaf + DMGreenLeaf;
TotDMDeadLeaf = TotDMDeadLeaf + DMDeadLeaf;
TotLeafArea = TotLeafArea + LeafArea;
TotDMLeafStem = TotDMLeaf + TotDMStem;

%total leaf area per m^2
LAI = TotLeafArea * PLM2;

%leaf area index
LeafCoverage = interp1(LeafCoverageConversion(:,2),...
                LeafCoverageConversion(:,1),DVS)/ 100;
LAICov = TotLeafArea*PLM2 / LeafCoverage;

% ROOT DEVELOPMENT
%Results of measurements done on roots cauliflower gamma dist
rootdensparam1 =interp1(RootTeSum, z1_param,TeSum); %Find current parameter1
rootdensparam2 =interp1(RootTeSum, z2_param,TeSum); %FInd current parameter2
Rdens = gampdf(1:TLayer:columndepth, rootdensparam1, rootdensparam2);
maxindex = find(Rdens == max(Rdens));
minindex = find(Rdens< 0.0001); %Calculate depths with low root densities
TotRootDepth = minindex(min(find((minindex>maxindex)==1))); %find the first depth after maxindex that has a low root densiti

rootdensparamxy =interp1(RootTeSum,x_y_param, TeSum);
cumulative_dist = normcdf(RootLateralDistance,0 , rootdensparamxy);
index = min(find(cumulative_dist > 1-Threshold));
fraction_plant = min(((RootLateralDistance(index)/100)^2*pi)/(1/PLM2)*1.2,1);
fraction_soil = 1- fraction_plant;

%Loop over layers to calculte sum of RDPF
NumLayer = round(TotRootDepth/TLayer);

RDMLA = zeros(85,1);
RLengthLA = zeros(85,1);

%Loop over layers
for i = 1:NumLayer
    RDMLA(i) = TotDMRoot * Rdens(i);
    RLengthLA(i) = RDMLA(i) * SRL;

end

% Calculate Radius of enclosing cilinder RCil
RCil = sqrt(fraction_plant/(PLM2*pi))*100;
RLDRCil = RLengthLA /((RCil^2) * pi * TLayer);

% CALCULATE NITROGEN DEMAND
NLeafStemOpt = interp1(CritN(:,2), CritN(:,1), DVS);
NShootOpt = (NLeafStemOpt * TotDMLeafStem + TotDMCurd * CurdNc + TotDMRoot * RootNc)/100;

demandoption = 1;
if demandoption ==1
    
    %N demand for optimal growth in leaves in g per plant
    LeafStemNDemand = NLeafStemOpt*...
        (TotDMGreenLeaf + TotDMStem)/100 - CumLeafStemNuptake;
    LeafStemNDemand = max(0, LeafStemNDemand);
    LeafNDemand = LeafStemNDemand*((TotDMGreenLeaf+TotDMDeadLeaf)/(TotDMGreenLeaf+TotDMDeadLeaf+TotDMStem));
    StemNDemand = LeafStemNDemand*((TotDMStem)/(TotDMGreenLeaf+TotDMDeadLeaf+TotDMStem));
    
    %N demand in the curds in g per plant
    CurdNDemand = CurdNc*(TotDMCurd/100) - CumCurdNuptake;
    CurdNDemand = max(CurdNDemand,0);
    
    %N demand in the roots in g per plant
    RootNDemand = RootNc*TotDMRoot/100-CumRootNuptake;
    RootNDemand = max(RootNDemand,0);

else
    %N demand for optimal growth in leaves in g per plant
    LeafStemNDemand = interp1(CritN(:,2),CritN(:,1),DVS) * (DMLeaf + DMStem) / 100;
    
    %N demand in the curds in g per plant
    CurdNDemand = CurdNc*DMCurd/100;
    
    %N demand in the roots in g per plant
    RootNDemand = RootNc*DMRoot/100;
end

%Calculate tota N demand per day in g per plant
ShootNDemand = LeafStemNDemand + CurdNDemand + RootNDemand;

%Integrate Tot N demand
TotNDemandCurd = TotNDemandCurd + CurdNDemand;
TotNDemandLeaf = TotNDemandLeaf + LeafNDemand;
TotNDemandLeafStem = TotNDemandLeafStem + LeafStemNDemand;
TotNDemandRoot = TotNDemandRoot + RootNDemand;
TotNDemandShoot = TotNDemandLeafStem + TotNDemandCurd + TotNDemandRoot;
TotNDemandStem = TotNDemandStem + StemNDemand;

% CALCULATE POTENTIAL TRANSPIRATION AND EVAPORATION ACCORDING TO FAO 56
idx = ET0_time==t;
epa = ET0_cm_per_day(idx)*(1-exp(-0.6*LAI));
esa = ET0_cm_per_day(idx)-epa;

CumEvaporation = CumEvaporation + esa;
CumTranspiration = CumTranspiration + epa;

%% FUNCTION OUTPUT

% crop state params
cropStateParams.JulianDay(end+1) = t; 
cropStateParams.CumCurdNuptake(end+1) = CumCurdNuptake; % cumulative curd nitrogen uptake [g/plant]
cropStateParams.CumEvaporation(end+1) = CumEvaporation; % cumulative evaporation [g/plant]
cropStateParams.CumLeafNuptake(end+1) = CumLeafNuptake; % cumulative leaf nitrogen uptake [g/plant]
cropStateParams.CumLeafStemNuptake(end+1) = CumLeafStemNuptake; % cumulative leaf and stem nitrogen uptake [g/plant]
cropStateParams.CumRootNuptake(end+1) =  CumRootNuptake; % cumulative root nitrogen uptake [g/plant]
cropStateParams.CumStemNuptake(end+1) =  CumStemNuptake; % cumulative stem nitrogen uptake [g/plang]
cropStateParams.CumShootNuptake(end+1) = CumShootNuptake; % cumulative shoot nitrogen uptake [g/plant]
cropStateParams.CumTranspiration(end+1) =  CumTranspiration; % cumulative transpiration
cropStateParams.CurdNDemand(end+1) =  CurdNDemand; % daily curd potential nitrogen demand [g/plant]
cropStateParams.DM_daily(end+1) = DM_daily;
cropStateParams.DMShoot(end+1) = DMShoot;% no param at the top
cropStateParams.DMCurd(end+1) =  DMCurd; % daily curd dry matter assimilation [g/plant]
cropStateParams.DMLeaf(end+1) =  DMLeaf; % daily leaf dry matter assimilation [g/plant]
cropStateParams.DMStem(end+1) =  DMStem; % daily stem dry matter assimilation [g/plant]
cropStateParams.DVS(end+1) =  DVS; % development stage
cropStateParams.LeafCoverage(end+1) =  LeafCoverage; % percentage soil covered by leaves
cropStateParams.LeafNDemand(end+1) =  LeafNDemand; % daily leaf nitrogen potential demand [g/plant]
cropStateParams.NShootOpt(end+1) = NShootOpt;
cropStateParams.Rdens(:,end+1) = Rdens;
cropStateParams.RootNDemand(end+1) =  RootNDemand; % daily root nitrogen potential demand [g/plant]
cropStateParams.ShootNDemand(end+1) =  ShootNDemand; % nitrogen demand of the whole plant [g/plant]
cropStateParams.StemNDemand(end+1) =  StemNDemand; % daily stem nitrogen potential demand [g/plant]
cropStateParams.TeSum(end+1) = TeSum; % cumulative effective temperature
cropStateParams.TotDM(end+1) = TotDM;
cropStateParams.TotDMCurd(end+1) = TotDMCurd; % total curd dry matter [g/plant]
cropStateParams.TotDMDeadLeaf(end+1) = TotDMDeadLeaf; % total dry matter of dead leaves [g/plant]
cropStateParams.TotDMGreenLeaf(end+1) = TotDMGreenLeaf; % total dry matter of green leaves [g/plant]
cropStateParams.TotDMLeaf(end+1) = TotDMLeaf; % total dry matter of leaves [g/plant]
cropStateParams.TotDMLeafStem(end+1) = TotDMLeafStem;
cropStateParams.TotDMRoot(end+1) = TotDMRoot; % total dry matter of the roots [g/plant]
cropStateParams.TotDMStem(end+1) = TotDMStem; % total dry matter of the stem [g/plant]
cropStateParams.TotNDemandCurd(end+1) = TotNDemandCurd; % total potential curd nitrogen demand [g/plant]
cropStateParams.TotNDemandLeaf(end+1) = TotNDemandLeaf;
cropStateParams.TotNDemandLeafStem(end+1) = TotNDemandLeafStem;
cropStateParams.TotNDemandRoot(end+1) = TotNDemandRoot; % total potential root nitrogen demand [g/plant]
cropStateParams.TotNDemandStem(end+1) = TotNDemandStem; % total potential stem nitrogen demand [g/plant]
cropStateParams.CurdNuptake(end+1) = CurdNuptake;
cropStateParams.DAP(end+1) = DAP;
cropStateParams.GP(end+1) =  GP; % gross photosynthetic production [g/plant]
cropStateParams.LeafStemNDemand(end+1) = LeafStemNDemand;
cropStateParams.NFac(end+1) = NFac;
cropStateParams.Novershoot(end+1) = Novershoot;
cropStateParams.RootNuptake(end+1) = RootNuptake;
cropStateParams.ShootNuptake(end+1) = ShootNuptake;
cropStateParams.StemNuptake(end+1) = StemNuptake;
cropStateParams.Te(end+1) = Te;
cropStateParams.TotLeafArea(end+1) = TotLeafArea;
cropStateParams.TotNDemandShoot(end+1) = TotNDemandShoot; %total potential shoot nitrogen demand [g/plant]
cropStateParams.TotRootDepth(end+1) = TotRootDepth; % total rooting depth [cm]
cropStateParams.fraction_plant(end+1) = fraction_plant;
cropStateParams.fraction_soil(end+1) = fraction_soil;
cropStateParams.LAI(end+1) = LAI; 
cropStateParams.LAICov(end+1) = LAICov; 
cropStateParams.LeafNuptake(end+1) = LeafNuptake;
cropStateParams.RCil(:,end+1) = RCil;
cropStateParams.RMAINT(end+1) = RMAINT; % daily maintenance respiration [g/plant]
cropStateParams.TotDMShoot(end+1) = TotDMShoot; % total dry matter of the whole plant [g/plant]
cropStateParams.RDMLA(:,end+1) =  RDMLA; % root dry matter per soil layer [g/plant]
cropStateParams.RLengthLA(:,end+1) = RLengthLA; % root length per soil layer [cm]
cropStateParams.RLDRCil(:,end+1) = RLDRCil;

% soil inner state params
soilInnerStateParams.epa = epa;
soilInnerStateParams.esa = 0;

% soil outer state params
soilOuterStateParams.epa = 0;
soilOuterStateParams.esa = esa;
