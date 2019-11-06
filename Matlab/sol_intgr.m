function [soilConsParams,soilStateCilinderParams] = sol_intgr(managementSettings,simulationSettings,...
    soilConsParams,soilStateCilinderParams)
%% DOCUMENTATION
%Calculates the mass balance for solute
%
%tsoli= inital mass solute in the profile
%tflsol total solute inflow in the soil surface
%rleasa total solute inflow at the bottom
%tsink= total sink term

%% FUNCTION INPUT

% managementSettings
nsol = managementSettings.nsol;
%simulation settings
dx = simulationSettings.dx;
immobile = simulationSettings.immobile;
ncomp = simulationSettings.nComp;

% soil constant parameters
alf = soilConsParams.alf;
bd = soilConsParams.bulkDensity;
tcsolo_ini = soilConsParams.tcsolo_ini;

acsolmo = soilStateCilinderParams.acsolmo;
acsolio = soilStateCilinderParams.acsolio;
acsolo = soilStateCilinderParams.acsolo;
cm = soilStateCilinderParams.cm;
cim = soilStateCilinderParams.cim;
cmp = soilStateCilinderParams.cmp;
cimp = soilStateCilinderParams.cimp;
csol = soilStateCilinderParams.csol;
dt = soilStateCilinderParams.dt;
initsol = soilStateCilinderParams.initsol;
rleasa = soilStateCilinderParams.rleasa;
sflxsa = soilStateCilinderParams.sflxsa;
sol_sinki = soilStateCilinderParams.sol_sinki;
sol_sinkm = soilStateCilinderParams.sol_sinkm;
tcsink = soilStateCilinderParams.tcsink;

tflsol = soilStateCilinderParams.tflsol;
tsoli = soilStateCilinderParams.tsoli;

wcmo = soilStateCilinderParams.wcmo;
wcio = soilStateCilinderParams.wcio;
wco = soilStateCilinderParams.wco;

%% FUNCTION MAIN BODY

volume=dx;
if initsol == 2
    joachim = 1;
    initsol = 0;
else
    joachim = 0;
end
if initsol==0
    tsoli=sum(tcsolo_ini);
    initsol = 1;
end 
ptsol=zeros(1,nsol);
ptsoli=zeros(1,nsol);
ptsolm=zeros(1,nsol);
ptsinki=zeros(1,nsol);
ptsinkm=zeros(1,nsol);
for j = 1:nsol  %%     do sp = 1,nr_of_sol
    %		calculate new total values of solute on the nodes
    for i=1:ncomp    %%%		do i=1,ncomp
        if immobile
            tcsolmo(i,j)=(cm(i,j)*wcmo(i)+(acsolmo(i,j)*...
                bd((i))))*volume;
            tcsolio(i,j)=(cim(i,j)*wcio(i)+(acsolio(i,j)*...
                bd((i))))*volume;
            tcsolo(i,j)=(tcsolmo(i,j)+tcsolio(i,j));
        else
            tcsolo(i,j)=(csol(i,j)*wcmo(i)+(acsolo(i,j)*...
                bd((i))))*volume;
            tcsolmo(i,j)=tcsolo(i,j);
            tcsolio(i,j)=0.0;
        end
    end
    
    
    for i=1:ncomp
        ptsoli(j)=ptsoli(j)+tcsolio(i,j);
        ptsolm(j)=ptsolm(j)+tcsolmo(i,j);
        ptsol(j)=ptsol(j)+tcsolo(i,j);
        tsinki(i,j) = (sol_sinki(i,j))*volume*dt;
        tsinkm(i,j) = (sol_sinkm(i,j)+alf((i))*...
            (cm(i,j)-cim(i,j)))*volume*dt;
        ptsinki(j) = ptsinki(j) +(sol_sinki(i,j))*volume*dt;
        
    end
    
    
    ptsinkm(j)= sum(sol_sinkm(:,j)+(alf.*((cm(1:ncomp,j)+cmp(1:ncomp,j))/2-(cim(:,j)+...
        cimp(:,j))/2)))*volume*dt;
    
    ptsink(j)=ptsinki(j)+ptsinkm(j) ;
    
    
    
    %      calculate cumulative values for the whole profile
    %		netto flux
    rleasa(j)=rleasa(j)+sflxsa(ncomp+1,j);
    %		cumulative downward flux (negativ)
    %if sflxsa(ncomp+1,j)<=0.0
    %dleasa(sp)=dleasa(j)+sflxsa(ncomp+1,j);
    %end
    %		cumulativ upward flux (positiv)
    %if sflxsa(ncomp+1,j)>=0.0
    %pleasa(j)=pleasa(j)+sflxsa(ncomp+1,j);
    %end
    %		inflow
    tflsol(j) = tflsol(j)+sflxsa(1,j);
    %		inflow during one day
    %	if(t- dint(t).eq.0.) solinfl(sp)=0.0
    %	solinfl(sp)=solinfl(sp)+sflxsa(1,sp)
    %		sink first order
    tcsink(j) = tcsink(j) + ptsink(j);
    %		change in the system from the start till the present time step
    dsol(j)=ptsol(j)-tsoli(j);
    %		mass balance error
    cberr_sol(j)=dsol(j) - (rleasa(j)+tflsol(j)+tcsink(j));
end

if  joachim == 1
    initsol =2;
end

%% FUNCTION OUTPUT

% soil constant parapemters
soilConsParams.tcsolo_ini = tcsolo_ini;

% soil stata cilinder parameters
soilStateCilinderParams.cberr_sol = cberr_sol;
soilStateCilinderParams.dsol = dsol;
soilStateCilinderParams.initsol = initsol;
soilStateCilinderParams.rleasa = rleasa;
soilStateCilinderParams.tcsink = tcsink;
soilStateCilinderParams.tcsolo = tcsolo;
soilStateCilinderParams.tflsol = tflsol;
soilStateCilinderParams.tsoli = tsoli;

        