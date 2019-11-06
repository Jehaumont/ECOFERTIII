function soilStateCilinderParams= solve_solute(climateState,cropConsParams,cropStateParams,managementSettings,...
    simulationSettings,soilConsParams,soilStateCilinderParams)
%% DOCUMENTATION
% Solve the 1-D convection dispersion equation
%using finite differences (cf. WAVE MANUEL)
%%%IN
%t= time(min)
%dt = the time increment (min)
%temp=temperature (°C)
%nsol=number of solute spieces
%WC:water content at time t, wco: water content=WC
%wat_flxs:(water fluxs ;ph (pressure head)
%solute_param1:some parameters of the solutes [bd,lm,ratio,alf,f]
%solute_param2:some parameters of the solutes [rates,kd,dif,ar,br]
%nitro_param:parameters for nitrogen simulation (from script in_nitro_om_param)
%soil_parameters: parameters of the soil
%om_param: parameters of the organic matter (from script in_nitro_om_param)
%rtex : root water extraction calculated in the water flow module
%wat_flxsa : water_flxs between the nodes (from script watsol)
%wat_flxsah:water_flxs between the nodes at half time (from script watsol)
%wcio: water content in the immobile region  (from script watso)
%wciob:water content in the immobile region at time t-dt (from script watsol)
%wcma:water content in the mobile region between nodes (from script watsol)
%wcmah: water content in the mobile region between nodes at half time (from script watsol)
%wcmo: water content in the mobile region at the nodes (from script watsol)
%wcmob:water content in the mobile region at the nodes at time t-dt (from script watsol)
%wco: water content at the nodes (from script watsol)
%pvela: pore velocity between the nodes (from script watsol)
%pvelah: pore velocity between the nodes at half time (from script watsol)
%pveloh: pore velocity at the nodes at half time (from script watsol)
%
%%%OUT
%cm:solute mobile concentration at time t
%cim:solute immobile concentration%at time t
%cs:solute concentration at the soil surface
%c_manp,n_manp,c_litp,n_litp,c_hump,n_hump: carbone and nitrogen in manure,
%litter and humus at time t
%diffus: diffusivity
%sol_sinkm: solute sink for the mobile water
%tot_upt: cumulative nitrogen uptake
%
%CALL:solute_boundary, solute_sink, miner_immob, nitro_uptake, nitro_sink
%CALLED BY wave_mat main programme
%-------------------------------------------------
% M. Sall 17/03/09

%% FUNCTION INPUT
% climate data
theta_table = climateState.climateDaily;

% simulation settings
dx = simulationSettings.dx;
immobile = simulationSettings.immobile;
ncomp = simulationSettings.nComp;
ncrop = simulationSettings.ncrop;
simplant = simulationSettings.simPlantFlag;
sim_nitro = simulationSettings.simNitroFlag;
t = simulationSettings.t;
tmax = simulationSettings.tmax;

% crop state parameters
drz = cropStateParams.TotRootDepth(end);
fraction_plant = cropStateParams.fraction_plant(end);
NDemand = cropStateParams.ShootNDemand(end);
rdens = cropStateParams.RLDRCil(:,end);
RLengthLA = cropStateParams.RLengthLA(:,end);

% crop constant parameters
plant_uptake_param = [cropConsParams.rorad,cropConsParams.rd0,cropConsParams.g,...
                     0,cropConsParams.rdens0,0];
                 
crop_type = cropConsParams.crop_type;
PropFact = [];

% management settings
harvest_date = sort([managementSettings.hDateCauli;managementSettings.hDateLeek]);
nsol = managementSettings.nsol;
plant_date = sort([managementSettings.pDateCauli;managementSettings.pDateLeek]);
if crop_type ==1
    PLM2 = managementSettings.PLM2Cauli;
elseif crop_type == 2
    PLM2 = managementSettings.PLM2Leek;
elseif isnan(crop_type)
    PLM2 = [];
end
ppdepth = managementSettings.ppdepth;
strategy = managementSettings.fertStrategy;

% soil constant parameters
bd = soilConsParams.bulkDensity;
alf = soilConsParams.alf;
ar = soilConsParams.ar;
br = soilConsParams.br;
ddepsol = soilConsParams.ddepsol;
dif = soilConsParams.dif;
f = soilConsParams.f;
kd = soilConsParams.kd;
lm = soilConsParams.lm;
lm(ncomp+1) = lm(ncomp);
miner_param = [soilConsParams.r0,soilConsParams.fe,soilConsParams.fh];
om_param = [soilConsParams.knitrif,soilConsParams.kdenit,soilConsParams.khydro,...
            soilConsParams.kvolat,soilConsParams.klit,soilConsParams.kman,soilConsParams.khum];
rates = soilConsParams.rates;
ratio = soilConsParams.ratio;
rcsolo = soilConsParams.rcsolo;
soil_parameters = [soilConsParams.wcr,soilConsParams.wcs,soilConsParams.alfa,...
                  soilConsParams.N,soilConsParams.ks,soilConsParams.lambda,soilConsParams.alfa_r];
solute_param1 = [soilConsParams.bulkDensity,soilConsParams.lm,soilConsParams.ratio,...
                soilConsParams.alf,soilConsParams.f];
solute_param2 = [soilConsParams.rates,soilConsParams.kd,soilConsParams.dif,soilConsParams.ar,...
                soilConsParams.br];          
tcsolo_ini = soilConsParams.tcsolo_ini;
wdepsol = soilConsParams.wdepsol;


% soil cilinder state params
acsolmo = soilStateCilinderParams.acsolmo;
acsolio = soilStateCilinderParams.acsolio;
cm = soilStateCilinderParams.cm;
conirsol = soilStateCilinderParams.conirsol;
cim = soilStateCilinderParams.cim;
csol = soilStateCilinderParams.csol;
cs = soilStateCilinderParams.cs;
diffush = soilStateCilinderParams.diffush;
dt = soilStateCilinderParams.dt;
first_time = soilStateCilinderParams.first_time;
fsol = soilStateCilinderParams.fsol;
initsol = soilStateCilinderParams.initsol;
om_appl = soilStateCilinderParams.om_appl;
ph = soilStateCilinderParams.ph;
pvela = soilStateCilinderParams.pvela;
pvelah = soilStateCilinderParams.pvelah;
pveloh = soilStateCilinderParams.pveloh;
reservoir = soilStateCilinderParams.reservoir;
rtex = soilStateCilinderParams.rtex;
soil_om = soilStateCilinderParams.soil_om;
tcsolo = soilStateCilinderParams.tcsolo;
temp = soilStateCilinderParams.temp;
tot_upt = soilStateCilinderParams.tot_upt;
uptakem = soilStateCilinderParams.uptakem;
uptakei = soilStateCilinderParams.uptakei;
wat_flxs = soilStateCilinderParams.wat_flxs;
wat_flxsa = soilStateCilinderParams.wat_flxsa;
wat_flxsah = soilStateCilinderParams.wat_flxsah;
WC = soilStateCilinderParams.WC;
wcio = soilStateCilinderParams.wcio;
wciob = soilStateCilinderParams.wciob;
wcma = soilStateCilinderParams.wcma;
wcmah = soilStateCilinderParams.wcmah;
wcmo = soilStateCilinderParams.wcmo;
wcmob = soilStateCilinderParams.wcmob;
wco = soilStateCilinderParams.wco;
wcob = soilStateCilinderParams.wcob;

%% FUNCTION MAIN BODY

%calculate beta coeficients
i = 1:ncomp;
betha1(i,1) = wat_flxsa(i)<0;
betha4(i,1) = 1 - betha1(i,1);
betha2(i,1)= wat_flxsa(i+1) <= 0 ;
betha3(i,1) = 1 - betha2(i,1);

%Determine initial values in the first time step
[cm,cim,csol,acsolmo,acsolio,initsol] = solute_data (cm,cim,csol,acsolmo,...
    acsolio,t,dx,tcsolo_ini,wco,wcmo,wcma,wcio,immobile,solute_param1,...
    solute_param2,nsol,ncomp,pvela,initsol);

%Determine the sink term
[decsolm,decsoli,sol_sinkm,sol_sinki] = solute_sink(dt,nsol,ncomp,cm,cim,wcmo,...
    wcio,solute_param1,solute_param2,immobile);


%%Nitrogen module

if sim_nitro
    %     [minerm,mineri,rcarbmin,soil_om]=...
    %      miner_immob(t,dt,temp,ncomp,WC,csol,ph,soil_om,wco,dx,cm,cim,soil_parameters,nitro_param,om_param);
    [minerm,mineri,rcarbmin,soil_om,om_appl,carbman, rnitman, carblit, rnitlit]=...
        miner_immob_new(t,dt,temp,ncomp,WC,csol,ph,soil_om,wco,dx,cm,cim,soil_parameters,om_param, miner_param, plant_uptake_param,om_appl,initsol);
    
    if simplant
        [uptakem,uptakei,tot_upt, first_time, uptake_matrix,rdens,und,unc]=nitro_uptake(csol,dt,dx,harvest_date,...
            ncomp,nsol,soil_parameters,plant_uptake_param,solute_param1,solute_param2,...
            plant_date,simplant,t,rtex,tot_upt,wco,first_time,drz,NDemand,PLM2,PropFact,...
            RLengthLA,rdens,ncrop,fraction_plant,rcsolo);
       
    else
        uptakem=zeros(ncomp,3);
        uptakei=zeros(ncomp,3);
        tot_upt=[];
        uptake_matrix = zeros(1,6);
        rdens= [];
        unc = zeros(ncomp,3);
        und = zeros(ncomp,3);
    end
    
    [sol_sink,sol_sinkm,sol_sinki,nitrifm,nitrifi,hydro_uream,hydro_ureai,uptakem,uptakei,denitm,...
        deniti,volatm,volati,decnorg,deccorg]=...
        nitro_sink(t,dt,tmax,temp,nsol,ncomp,cm,cim,WC,wcmo,wcio,decsolm,decsoli,minerm,mineri,rcarbmin,...
        uptakem,uptakei,soil_parameters,solute_param1,solute_param2,om_param, miner_param, plant_uptake_param,csol,immobile);
    
else
    uptake_matrix = zeros(1,6);
    minerm= [];
    mineri=[];
    nitrifm=[];
    nitrifi=[];
    denitm=[];
    deniti=[];
    hydro_uream=[];
    hydro_ureai=[];
    volatm=[];
    volati=[];
    decnorg=[];
    deccorg=[];
    carbman=[];
    carblit=[];
    rnitman=[];
    rnitlit=[];
    cco2o=[];
    drza=[];
end

% actual value to passed value
[conc,solute_applic,solboco_top_type,solboco_input_type]=In_solute_boundary(t,nsol);
if solboco_input_type ==2
    csp=cs./abs(wat_flxsah(1))/dt;
    if sum(isnan(csp))== nsol
        csp = zeros(1,nsol);
    end
else
    csp =cs;
end
cmp=cm;
cimp=cim; % cm and cimp = precedent mobile and immobile soute concentration
acsolmb=acsolmo;
%set cmp zero gradient
cmp(ncomp+1,:)=cmp(ncomp,:);

% Boundary Condition
%[conc,solute_applic,solboco_top_type,solboco_input_type]=In_solute_boundary(t,nsol);
if solboco_input_type == 1 % As by M. Sall expressed in M L-3
    if wat_flxs(1)<=0
        cs=conc;
    else
        cs=zeros(1,nsol);
    end
elseif solboco_input_type ==2 %As by J. Vansteenkiste expressed in M L-2
    
    if wat_flxs(1) <0
        [fsol, conirsol, solsur] = solute_upper(t, fsol, conirsol, wdepsol, ddepsol, dt, strategy, theta_table, fraction_plant);
        cs = (reservoir + solsur)./abs(wat_flxsah(1))/dt;
        cs_temp = (reservoir + solsur);
        reservoir = zeros(1,nsol);
        if sum(cs) ~=0
            disp('')
        end
    else
        solsur =[0 0 0];
        reservoir = (reservoir + solsur);
        if sum(reservoir) ~=0
            %disp('')
        end
        cs =zeros(1,nsol);
        cs_temp =cs;
    end
end

%Determine solute redistribution due to plowing
if solboco_input_type == 2
    [csol, cm, cim, acsolmo, acsolio,ppdepth,tcsolo] = plowing(t,solute_param1, solute_param2,csol, cm, cim, acsolmo,...
        acsolio, wcio, wco, nsol, immobile,wcmo,tcsolo,dx,ppdepth,solsur);
end


cm = inv_sol_matrix(immobile,ncomp,cs,f,dt, pveloh, wat_flxsah, wcmo, wcmob, dx, diffush, wcmah, pvelah,betha1,betha2, betha3,betha4, csp, cmp, bd, kd, sol_sinkm,nsol);

for j=1:nsol
    %%%%
    for I=1:ncomp
        if cm(I,j)<0
            cm(I,j)=0;
        end
        %Adsorbed mass
        acsolmo(I,j)=kd(I,j)*f(I)*cm(I,j);
        %SOLUTE IMMOBILE CONCENTRATION
        if immobile
            ag(I,j)=wciob(I)/dt+(1-f(I)).*bd(I).*kd(I,j)/dt-alf(I)/2;
            bg(I,j)=wcio(I)/dt+(1-f(I))*bd(I)*kd(I,j)/dt+alf(I)/2;
            cim(I,j)=(ag(I,j)*cimp(I,j)+(alf((I))*cm(I,j))/2 +...
                (alf((I))*cmp(I,j))/2)/bg(I,j)+sol_sinki(I,j)/bg(I,j);
            acsolio(I,j)=cim(I,j)*kd((I),j)*(1-f((I)));
            csol(I,j)=(cim(I,j)*wciob(I)+cm(I,j)*wcmob(I))/wcob(I);
            acsolo(I,j)=acsolio(I,j)+acsolmo(I,j);
        else
            cim(I,j)=0.0;
            acsolio(I,j)=0.0;
            csol(I,j)=cm(I,j);
            acsolo(I,j)=acsolmo(I,j);
        end
    end
    
    for i=2:ncomp
        if wat_flxsah(i)<0
            sflxsa(i,j)=(wat_flxsah(i)*(cm((i-1),j)+...
                cmp((i-1),j))/2)*dt;
        else
            sflxsa(i,j)=(wat_flxsah(i)*(cm(i,j)+...
                cmp(i,j))/2)*dt;
        end
    end
    if wat_flxsah(1)<0
        sflxsa(1,j) = (cs(j)+csp(j))/2*abs(wat_flxsah(1))*dt+(diffush(1,j))*...
            ((cs(j)+csp(j))-(cm(1,j)+cmp(1,j)))/2/dx*(wcmah(1))*dt;
    else
        sflxsa(1,j)=0.0;
    end
    volume= dx;
    if  abs(cm(ncomp,j)-cmp(ncomp,j)>0.001  & wat_flxsah(ncomp+1)>0.00000000001)
        sflxsa((ncomp+1),j)=((-cmp(ncomp,j)+cm(ncomp,j))*...
            wcmo(ncomp)*dx)-(acsolmb(ncomp,j)-acsolmo(ncomp,j))*...
            bd(1)*volume;
    elseif wat_flxsah(ncomp+1)<0.0
        sflxsa((ncomp+1),j)=(wat_flxsah(ncomp+1)*...
            ((cm(ncomp,j)+ cmp(ncomp,j))/2))*dt;
    else
        sflxsa(ncomp+1,j)=(wat_flxsah(ncomp+1)*cmp((ncomp+1),j))*dt;
    end
end

if solboco_input_type ==2
    cs= cs_temp;
end

%% FUNCTION OUTPUT

% soil state cilinder params
soilStateCilinderParams.acsolio = acsolio;
soilStateCilinderParams.acsolmo = acsolmo;
soilStateCilinderParams.acsolo = acsolo;
soilStateCilinderParams.carblit = carblit;
soilStateCilinderParams.carbman = carbman;
soilStateCilinderParams.cim = cim;
soilStateCilinderParams.cimp = cimp;
soilStateCilinderParams.cm = cm;
soilStateCilinderParams.cmp = cmp;
soilStateCilinderParams.conirsol = conirsol;
soilStateCilinderParams.cs = cs;
soilStateCilinderParams.csol = csol;
soilStateCilinderParams.csp = csp;
soilStateCilinderParams.ddepsol = ddepsol;
soilStateCilinderParams.deccorg = deccorg;
soilStateCilinderParams.decnorg = decnorg;
soilStateCilinderParams.deniti = deniti;
soilStateCilinderParams.denitm = denitm;
soilStateCilinderParams.first_time = first_time;
soilStateCilinderParams.fsol = fsol;
soilStateCilinderParams.hydro_ureai = hydro_ureai;
soilStateCilinderParams.hydro_uream = hydro_uream;
soilStateCilinderParams.mineri = mineri;
soilStateCilinderParams.minerm = minerm;
soilStateCilinderParams.nitrifi = nitrifi;
soilStateCilinderParams.nitrifm = nitrifm;
soilStateCilinderParams.om_appl = om_appl;
soilStateCilinderParams.ppdepth = ppdepth;
soilStateCilinderParams.rdens = rdens;
soilStateCilinderParams.reservoir = reservoir;
soilStateCilinderParams.rnitlit = rnitlit;
soilStateCilinderParams.rnitman = rnitman;
soilStateCilinderParams.sflxsa = sflxsa;
soilStateCilinderParams.soil_om = soil_om;
soilStateCilinderParams.sol_sinki = sol_sinki;
soilStateCilinderParams.sol_sinkm = sol_sinkm;
soilStateCilinderParams.solute_applic = solute_applic;
soilStateCilinderParams.tot_upt = tot_upt;
soilStateCilinderParams.unc = unc;
soilStateCilinderParams.und = und;
soilStateCilinderParams.uptake_matrix = uptake_matrix;
soilStateCilinderParams.uptakei = uptakei;
soilStateCilinderParams.uptakem = uptakem;
soilStateCilinderParams.volati = volati;
soilStateCilinderParams.volatm = volatm;
soilStateCilinderParams.wdepsol = wdepsol;
