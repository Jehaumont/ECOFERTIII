function [cropStateParams,soilInnerStateParams,soilOuterStateParams] =...
                                            run_model_leek(climateConsParams,...
                                            climateState,cropStateParams,leekConsParams,...
                                            fileSettings, managementSettings,...
                                            simulationSettings,soilCommonStateParams,...
                                            soilInnerStateParams,soilOuterStateParams)
%% FUNCTION INPUT
% climate constant parameters
AtmP = climateConsParams.atmp; % atmospheric pressure 
Cp = climateConsParams.Cp; %specific heat of air
CO2L = climateConsParams.CO2L; %CO2 level

% climate data
ET0_time = climateState.ET0_time;
ET0_cm_per_day = climateState.ET0_cm_per_day;
RefClim = climateState.crop_climate;


% management settings
ponsDepth = managementSettings.ponsDepthLeek;
if simulationSettings.calibrateFlag ==1
    PLM2 = 18.18;
else
    PLM2 = managementSettings.PLM2Leek; % plant density [plants/M?]
end

plant_date = sort([managementSettings.pDateCauli;managementSettings.pDateLeek]);

% soil common parameters 
NSupply = soilCommonStateParams.NSupply; % the amount of nitrogen [g] that is available for the plant to take up

% simulation settings
columndepth = simulationSettings.columnDepth; % depth of soil column for simulation
ncrop = simulationSettings.ncrop; % number of rotation that is currently on the field
NFAST = simulationSettings.Nfast; % iterations in the fast loop per day
t = simulationSettings.t; % time of simulation
TLayer = simulationSettings.dx; % thickness of a soil layer [cm]
DAP = t-plant_date(ncrop)+1;
DTFAST = 1.0/NFAST;

% leek constant parameters
albedo = leekConsParams.albedo;
CK = leekConsParams.CK; % exponential decay factor for vapour pressure deficit
DVSeffect = leekConsParams.DVSEffect; % reduction factor on maximal photosynthesis based on DVS [reduction factor,DVS]
Ilsolrad = leekConsParams.Ilsolrad; % conversion factor of radiation to long wave radiation
Issolrad = leekConsParams.Issolrad; % conversion factor of radiatino to short wave radiation
GREF = leekConsParams.GREF; % growth efficiency g dry matter per g assimilate
LeafCoverageConversion = leekConsParams.LeafCoverageConversion; % %soil coverd by leaves based on DVS [%coverage,DVS]
LeafFrac = leekConsParams.LeafFrac; % fractional allocation of assimilated dry matter to the leaves related to DVS [%,DVS]
PARsolrad = leekConsParams.PARsolrad; % photosynthetically active radiation
PPFDpar = leekConsParams.PPFDpar;
PercNLeaf = leekConsParams.PercNLeaf; % critical nitrogen concentration in the leaves
PercNShaft = leekConsParams.PercNShaft; % critical nitrogen concentration in the shaft
Q10 = leekConsParams.Q10; % increase of maintenance respiration when temperature increases 10?C
Qe = leekConsParams.Qe; % PAR utilisation coefficient [mgCO2/J]
RDPF = leekConsParams.RDPF; % root density distribution in depth axis
RGRDi = leekConsParams.RGRDi; % exponential coefficient for the conversion of DM to shaft diameter
RGRLe = leekConsParams.RGRLe; % exponential coefficient for the converison of DM to shaft length
RMRL = leekConsParams.RMRL; % respiration coefficient for the leaves [gCH20/(gDM*Day)]
RMRF = leekConsParams.RMRF; % respiration coefficient for the curd/stem [gCH20/(gDM*Day)]
RootFrac = leekConsParams.RootFrac; % percentage of dry matter allocated to the roots 
RootTeSum = leekConsParams.RootTeSum; % effective temperature sum for root development
rr = leekConsParams.rr; % root action radius [cm]
SLAEffect = leekConsParams.SLAEffect; % sla converstion to leaf area based on DVS [converstion factor, DVS]
SRL = leekConsParams.SRL; % specific root length [cm/g]
StemFrac = leekConsParams.StemFrac; % fractional allocation of assimilated dry matter to the stem related to DVS [%, DVS]
TAU = leekConsParams.TAU; % leaf conductance for C02 transfer [m/s]
Tb = leekConsParams.Tb; % base temperature [?C]
TempEffect = leekConsParams.TempEffect; % temperature reduction on maximal photosynthetic efficiency [reduction factor, temperature]
virRd = leekConsParams.virRd;
VPDL = leekConsParams.VPDL; % vapour pressure deficit limit [kPa]
XK = leekConsParams.XKK; % light extinsion coefficient
XM = leekConsParams.XM; % leaf transmission coefficient
x_y_param = leekConsParams.x_y_param; % parameter for the lateral distribution of the roots
z1_param = leekConsParams.z1_param; % gamma distribution parameter for estimation of root denstiy
z2_param = leekConsParams.z2_param; % gamma distribution parameter for estimation of root density

% leek state parameters
CumEvaporation = cropStateParams.CumEvaporation(end); % cumulative evaporation [g/plant]
CumLeafNuptake = cropStateParams.CumLeafNuptake(end); % cumulative leaf nitrogen uptake [g/plant]
CumLeafStemNuptake = cropStateParams.CumLeafStemNuptake(end);
CumRootNuptake = cropStateParams.CumRootNuptake(end); % cumulative root nitrogen uptake [g/plant]
CumStemNuptake = cropStateParams.CumStemNuptake(end); % cumulative stem nitrogen uptake [g/plang]
CumShootNuptake = cropStateParams.CumShootNuptake(end);% CumRootNuptake + CumStemNuptake + CumLeafNuptake; % cumulative shoot nitrogen uptake [g/plant]
CumTranspiration = cropStateParams.CumTranspiration(end); % cumulative transpiration
DMStem = cropStateParams.DMStem(end); % daily stem dry matter assimilation [g/plant]
LeafCoverage = cropStateParams.LeafCoverage(end); % percentage soil covered by leaves
LeafNDemand = cropStateParams.LeafNDemand(end); % daily leaf nitrogen demand [g/plant]
RootNDemand = cropStateParams.RootNDemand(end); % daily root nitrogen demand [g/plant]
ShootNDemand = cropStateParams.ShootNDemand(end); % daily nitrogen demand of the whole plant [g/plant]
StemNDemand = cropStateParams.StemNDemand(end); % daily stem nitrogen demand [g/plant]
TeSum = cropStateParams.TeSum(end); % cumulative effective temperature
TotDM = cropStateParams.TotDM(end); % total dry matter produced [g/plant]
TotDMLeaf = cropStateParams.TotDMLeaf(end); % total dry matter of leaves [g/plant]
TotDMLeafStem = cropStateParams.TotDMLeafStem(end); % total dry matter of leaf and stem [g/plant]
TotDMShoot = cropStateParams.TotDMShoot(end); % total dry matter of the whole plant [g/plant]
TotDMStem = cropStateParams.TotDMStem(end); % total dry matter of the stem [g/plant]
TotDMRoot = cropStateParams.TotDMRoot(end); % total dry matter of the roots [g/plant]
TotNDemandLeaf = cropStateParams.TotNDemandLeaf(end); % total leaf nitrogen demand
TotNDemandRoot = cropStateParams.TotNDemandRoot(end); % total root nitrogen demand [g/plant]
TotNDemandStem = cropStateParams.TotNDemandStem(end); % total stem nitrogen demand [g/plant]
TotNDemandShoot = cropStateParams.TotNDemandShoot(end);
TotNOptLeaf = cropStateParams.TotNOptLeaf(end);
TotNOptStem = cropStateParams.TotNOptStem(end);
DMLeaf = cropStateParams.DMLeaf(end); % no relocation, daily leaf dry matter assimilation [g/plant]

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
else
    ratioleaf = 0;
    ratiostem = 0;
    ratioroot = 0;
end

LeafNuptake = ShootNuptake*ratioleaf;
StemNuptake = ShootNuptake*ratiostem;
LeafStemNuptake = LeafNuptake + StemNuptake; 
RootNuptake = ShootNuptake*ratioroot;

% integration of total nitrogen uptake and calculation of Nfac
if DAP ~= 1
    
    % Integration of total nitrogen uptake
    CumLeafNuptake = CumLeafNuptake + LeafNuptake;
    CumStemNuptake = CumStemNuptake + StemNuptake;
    CumLeafStemNuptake = CumLeafStemNuptake + LeafStemNuptake;
    CumRootNuptake = CumRootNuptake + RootNuptake;
    CumShootNuptake = CumShootNuptake + ShootNuptake;

    % Calculation of nitrogen correction factor
    PercNLeaftab = CumLeafNuptake/TotDMLeaf * 100;
    PercNStemtab = CumStemNuptake/TotDMStem*100;
    
    NFac_Leaf = PercNLeaftab/interp1(PercNLeaf(:,2), PercNLeaf(:,1), TeSum);
    NFac_Stem = PercNStemtab/interp1(PercNShaft(:,2), PercNShaft(:,1), TeSum);    
    NFac_Leaf = min(1, NFac_Leaf);
    NFac_Stem = min(1, NFac_Stem);
    
    DMratio_Leaf = TotDMLeaf/TotDMLeafStem;
    DMratio_Stem = TotDMStem/TotDMLeafStem;
    
    % calculated weighted mean of both correction factors
    NFac = NFac_Leaf * DMratio_Leaf + NFac_Stem * DMratio_Stem;    
    NFac = min(NFac,1);
    NFac = max(NFac, 0.05);
end

% Initialize total transpiration and evaporation and aily temperature sum
GP =0;
RMAINT = 0;
Tsum = 0;

for JF = 1 : NFAST
    
    % TFAST is the hour of the day (h)
    TFAST = (JF-1)* 24/NFAST;
    
    d = find(RefClim(:,1)==t & RefClim(:,2)==(TFAST+1));
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
    Ps = 0.61078 * exp((17.2694*TMPA)/(TMPA+237.3));
    %Calculate actual vapor pressure PV
    Pv = RH/100 * Ps;
    VPD = Ps-Pv;
    PAR = solrad * PARsolrad;
    PPFD = PAR *PPFDpar;
    
    % CALCULATION OF THE GROSS PHOTOSYNTHETIC PRODUCTION
    GPF = 0;
    
    %Effect of CO2 on Pmax (from Gainesville)
    PMAX = TAU(1) * CO2L;
    if CO2L>1500
        PMAX = (TAU(1) * 1500) + (TAU(2) * (CO2L-1500));
    end
    %Reduction of Pmax at extreme temperatures
    PMAX = PMAX * interp1(TempEffect(:,2),TempEffect(:,1),TMPA);
    
    %Reduction of Pmax in relation to DVS, senescence  in LEEK
    PMAX = PMAX * interp1(DVSeffect(:,2),DVSeffect(:,1),TeSum);
    
    % Reduction of Pmax related to N stress
    PMAX = PMAX * NFac;
        
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
        % conversion of GPF from ?M/m2-s into g/m2-day
        %?M/m2-s x 0.000044g/?M x 3600s/h x 24h/d = 3.8016 g/m2-day
        % GPF in g per day per plant
        GPF = GPF * 3.8016/PLM2;
    end
    
    % MAINTENANCE RESPIRATION
    % Effect of temperature on maintenance respiration on hourly basis
    RMAINTF = (RMRL * (DMLeaf)+RMRF * DMStem) * (Q10^(0.1 * TMPA-2.0));
    
    %integration of variables on 24 hours
    GP = GP + GPF * DTFAST;
    RMAINT = RMAINT + RMAINTF * DTFAST;
        
    % CALUCULATE TEMPERATURE SUM
    Tsum = Tsum + TMPA;
    
end

% DAILY LOOP
%GP and RMAINT in g CH2O per plant
%Transfert to Biomass, calculate Total DM in g per plant per dag
DM_daily = (GP - RMAINT) * GREF; 

%Calculate T Sum  as measure for development stage
MeanDayTemp = Tsum/24;
Td = MeanDayTemp - Tb;
%Effective temperature set to zero if negative, non-decreasing Tsum
if Td > 0
    Te = Td;
else
    Te = 0;
end
%Calculate Temperature sum
TeSum = TeSum + Te;

%Compartimentation over leaf, stem and curd
%Growth rates g per planr per day
DMLeaf = DM_daily * interp1(LeafFrac(:,2),LeafFrac(:,1),max(TeSum, LeafFrac(end,2)))/(1 + RootFrac);
DMStem = DM_daily * interp1(StemFrac(:,2),StemFrac(:,1),max(TeSum, StemFrac(end,2)))/(1 + RootFrac);

%Compartimentation to roots, 9 % of shoot and constant
DMShoot = DMLeaf + DMStem;
DMRoot = DMShoot * RootFrac;

%Calculate leaf area out of SLA-DVS curve
LeafArea = DMLeaf * (interp1(SLAEffect(:,2),SLAEffect(:,1),TeSum)/10000);

%Integrate to cumulative total DM
TotDM =  TotDM + DM_daily;
TotDMLeaf = TotDMLeaf + DMLeaf;
TotDMStem = TotDMStem + DMStem;
TotDMLeafStem = TotDMLeaf + TotDMStem;
TotDMShoot = TotDMShoot + DMShoot;
TotDMRoot = TotDMRoot + DMRoot;
TotLeafArea = TotLeafArea + LeafArea;

%Diameter (mm) en length (cm) calculations
ShaftDi = exp(RGRDi(2)) * exp(RGRDi(1) * TeSum);
ShaftLe = exp(RGRLe(2)) * exp(RGRLe(1) * TeSum);

%Direct LAI calculation 
LAI = TotLeafArea * PLM2;
%direct LAI calculation per m^2 coverage
LeafCoverage = interp1(LeafCoverageConversion(:,2),LeafCoverageConversion(:,1),TeSum)/100;
LAICov = TotLeafArea * PLM2 / LeafCoverage;

% ROOT DEVELOPMENT
% soil root coverage
leek_rootcoverage = interp1(RootTeSum,x_y_param,TeSum);
leek_rootcoverage = max(leek_rootcoverage*1.4, 1);
fraction_plant = leek_rootcoverage/100;
fraction_soil = 1-fraction_plant;


%RootDepth in cm follows  linear curve in function of TeSum
rootdensparam1 = interp1(RootTeSum,z1_param,TeSum);
rootdensparam2 = interp1(RootTeSum,z2_param,TeSum);
Rdens = zeros(columndepth,1);
Rdens(ponsDepth+1:end) = gampdf(1:TLayer:columndepth-ponsDepth,rootdensparam1,rootdensparam2);
maxIndex = find(Rdens == max(Rdens));
minIndex = find(Rdens < 0.001);
TotRootDepth = minIndex(min(find((minIndex>maxIndex)==1)));

%Loop over layers to calculte sum of RDPF
NumLayer = round(TotRootDepth/TLayer);
RDMLA = zeros(85,1);
RLengthLA = zeros(85,1);

%Loop over layers
for i = 1 : NumLayer

    RDMLA(i) = TotDMRoot * Rdens(i);
    RLengthLA(i) = RDMLA(i) * SRL;

end

%Calculate Radius of enclosing cilinder RCil
%Root Length Density per layer in cm/cm^3 for enclosing cilinder
RCil = sqrt(fraction_plant/(PLM2*pi))*100;
RLDRCil = RLengthLA /((RCil^2) * pi * TLayer);

% CALCULATE NITROGEN UPTAKE FOR MAXIMAL GROWTH
%N in leaves
NLeafOpt = interp1(PercNLeaf(:,2),PercNLeaf(:,1),TeSum);
NStemOpt = interp1(PercNShaft(:,2),PercNShaft(:,1),TeSum);

%Integrate Tot N
TotNOptLeaf = NLeafOpt/100 * TotDMLeaf;
TotNOptStem = NStemOpt/100 * TotDMStem;
TotNOptLeafStem = TotNOptLeaf + TotNOptStem;

demandoption = 1;
if demandoption ==1
    %N demand for optimal growth in leaves in g per plant
    LeafNDemand = NLeafOpt * (TotDMLeaf)/100 - CumLeafNuptake;
    LeafNDemand = max(0,LeafNDemand);
    
    StemNDemand = NStemOpt * (TotDMStem)/100 - CumStemNuptake;
    StemNDemand = max(0,StemNDemand);
    LeafStemNDemand = StemNDemand + LeafNDemand;
    RootNDemand = 0;

else
    LeafNDemand = NperLeaf;
    StemNDemand = NperStem ;
    LeafStemNDemand = StemNDemand + LeafNDemand;
    RootNDemand = 0;
end

%Calculate tota N demand per day in g per plant
ShootNDemand = LeafNDemand + StemNDemand + RootNDemand;

%Integrate Tot N demand  
TotNDemandLeaf = TotNDemandLeaf + LeafNDemand;
TotNDemandStem = TotNDemandStem + StemNDemand;
TotNDemandRoot = TotNDemandRoot + RootNDemand;
TotNDemandShoot = TotNDemandShoot + ShootNDemand;

% CALCULATE POTENTIAL TRANSPIRATION AND EVAPORATION ACCORDING TO FAO 56
idx = ET0_time==t;
epa = ET0_cm_per_day(idx)*(1-exp(-0.6*LAI));
esa = ET0_cm_per_day(idx)-epa;

CumEvaporation = CumEvaporation + esa;
CumTranspiration = CumTranspiration + epa;

%% FUNCTION OUTPUT

cropStateParams.JulianDay(end+1) = t;
cropStateParams.CumEvaporation(end+1) = CumEvaporation; % cumulative evaporation [g/plant]
cropStateParams.CumLeafNuptake(end+1) = CumLeafNuptake; % yes top, cumulative leaf nitrogen uptake [g/plant]
cropStateParams.CumLeafStemNuptake(end+1) = CumLeafStemNuptake; % yes top
cropStateParams.CumRootNuptake(end+1) = CumRootNuptake; % yes top, cumulative root nitrogen uptake [g/plant]
cropStateParams.CumStemNuptake(end+1) = CumStemNuptake; % yes top, cumulative stem nitrogen uptake [g/plang]
cropStateParams.CumShootNuptake(end+1) = CumShootNuptake; % yes top, cumulative shoot nitrogen uptake [g/plant]
cropStateParams.CumTranspiration(end+1) = CumTranspiration; % yes top,cumulative transpiration
cropStateParams.DM_daily(end+1) = DM_daily;
cropStateParams.DMLeaf(end+1) = DMLeaf; % no top daily leaf dry matter assimilation [g/plant]
cropStateParams.DMShoot(end+1) = DMShoot; %no top
cropStateParams.DMStem(end+1) = DMStem; % yes top, daily stem dry matter assimilation [g/plant]
cropStateParams.LeafCoverage(end+1) = LeafCoverage; % yes top, percentage soil covered by leaves
cropStateParams.LeafNDemand(end+1) = LeafNDemand; % yes top, daily leaf nitrogen demand [g/plant]
cropStateParams.RootNDemand(end+1) = RootNDemand; % yes top, daily root nitrogen demand [g/plant]
cropStateParams.ShootNDemand(end+1) = ShootNDemand; % yes top, daily nitrogen demand of the whole plant [g/plant]
cropStateParams.StemNDemand(end+1) = StemNDemand; % yes top, daily stem nitrogen demand [g/plant]
cropStateParams.TeSum(end+1) = TeSum; % yes top, cumulative effective temperature
cropStateParams.TotDM(end+1) = TotDM;
cropStateParams.TotDMLeaf(end+1) = TotDMLeaf; % yes top, total dry matter of leaves [g/plant]
cropStateParams.TotDMLeafStem(end+1) = TotDMLeafStem;
cropStateParams.TotDMShoot(end+1) = TotDMShoot; % yes top, total dry matter of the whole plant [g/plant]
cropStateParams.TotDMStem(end+1) = TotDMStem; % yes top, total dry matter of the stem [g/plant]
cropStateParams.TotDMRoot(end+1) = TotDMRoot; % yes top, total dry matter of the roots [g/plant]
cropStateParams.TotNDemandLeaf(end+1) = TotNDemandLeaf; % yes top, total leaf nitrogen demand
cropStateParams.TotNDemandRoot(end+1) = TotNDemandRoot; % yes top, total root nitrogen demand [g/plant]
cropStateParams.TotNDemandStem(end+1) = TotNDemandStem; % yes top, total stem nitrogen demand [g/plant]
cropStateParams.TotNDemandShoot(end+1) = TotNDemandShoot; % yes top,
cropStateParams.TotNOptLeaf(end+1) = TotNOptLeaf; % yes top
cropStateParams.TotNOptStem(end+1) = TotNOptStem; % yes top
cropStateParams.TotNOptLeafStem(end+1) = TotNOptLeafStem;
cropStateParams.DAP(end+1) = DAP; %no top
cropStateParams.RDMLA(:,end+1) = RDMLA; % no top, root dry matter per soil layer [g/plant] 
cropStateParams.RLengthLA(:,end+1) = RLengthLA; % no top, root length per soil layer [cm]
cropStateParams.Rdens(:,end+1) = Rdens;
cropStateParams.RLDRCil(:,end+1) = RLDRCil;
cropStateParams.GP(end+1) = GP; % gross photosynthetic production [g/plant]
cropStateParams.LAI(end+1) = LAI;
cropStateParams.LAICov(end+1) = LAICov;
cropStateParams.LeafNuptake(end+1) = LeafNuptake;
cropStateParams.LeafStemNDemand(end+1) = LeafStemNDemand;
cropStateParams.NFac(end+1) = NFac;
cropStateParams.Novershoot(end+1) = Novershoot;
cropStateParams.RCil(end+1) = RCil;
cropStateParams.RMAINT(end+1) = RMAINT; % daily maintenance respiration [g/plant]
cropStateParams.RootNuptake(end+1) = RootNuptake;
cropStateParams.ShaftDi(end+1) = ShaftDi;
cropStateParams.ShaftLe(end+1) = ShaftLe;
cropStateParams.ShootNuptake(end+1) = ShootNuptake;
cropStateParams.StemNuptake(end+1) = StemNuptake;
cropStateParams.Te(end+1) = Te;
cropStateParams.TotLeafArea(end+1) = TotLeafArea;
cropStateParams.TotRootDepth(end+1) = TotRootDepth;
cropStateParams.fraction_plant(end+1) = fraction_plant;
cropStateParams.fraction_soil(end+1) = fraction_soil;

% soil inner state params
soilInnerStateParams.epa = epa;
soilInnerStateParams.esa = 0;

% soil outer state params
soilOuterStateParams.epa = 0;
soilOuterStateParams.esa = esa;
