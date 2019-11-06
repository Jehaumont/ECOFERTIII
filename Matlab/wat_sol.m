function soilStateCilinderParams = wat_sol(managementSettings,simulationSettings,soilConsParams,soilStateCilinderParams)
%% DOCUMENTATION
%Prepare data from the water flow module  for the solute module     
%mainly inspired from WATSOL.FOR         
%IN
%t:time
%wat_flxs:water fluxs 
%immobile: choice for immobile water concept (see in_solute_parameter)
%solute_param1: parameters of the solutes
%solute_param2: parameters of the solutes
%ncomp: number of compartiments
%WC: water content
%nsol:number of solutes
%
%OUT
%wat_flxsa : water_flxs between the nodes
%wat_flxsah:water_flxs between the nodes at half time
%wcio: water content in the immobile region 
%wciob:water content in the immobile region at time t-dt
%wcma:water content in the mobile region between nodes 
%wcmah: water content in the mobile region between nodes at half time
%wcmo: water content in the mobile region at the nodes 
%wcmob:water content in the mobile region at the nodes at time t-dt
%wco: water content at the nodes
%pvela: pore velocity between the nodes
%pvelah: pore velocity between the nodes at half time
%pveloh: pore velocity at the nodes at half time
%pvelo: pore velocity between the nodes 
%diffush: apparent diffusivity
%CALL: none
%CALLED BY: solve_solute
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% M. Sall 10/04/09

%% FUNCTION INPUT

% management settings
nsol = managementSettings.nsol;

% simulation settings
immobile = simulationSettings.immobile;
ncomp = simulationSettings.nComp;
t = simulationSettings.t;

% soil constant parameters
ar = soilConsParams.ar;
br = soilConsParams.br;
dif = soilConsParams.dif;
lm = soilConsParams.lm;
lm(ncomp+1) = lm(ncomp);
ratio = soilConsParams.ratio;

% soil state cilinder params
diffus = soilStateCilinderParams.diffus;
initsol = soilStateCilinderParams.initsol;
pvela = soilStateCilinderParams.pvela;
pvelo = soilStateCilinderParams.pvelo;
wat_flxs = soilStateCilinderParams.wat_flxs;
wat_flxsa = soilStateCilinderParams.wat_flxsa;
WC = soilStateCilinderParams.WC;
wcio = soilStateCilinderParams.wcio;
wcma = soilStateCilinderParams.wcma;
wcmo = soilStateCilinderParams.wcmo;
wco = soilStateCilinderParams.wco;

%% FUNCTION MAIN BODY

if initsol == 2
    joachim = 1;
    initsol = 0;
else
    joachim = 0;
end
if initsol == 0 %first time in this program
%if t==0 

%		calculation of the fluxes and water contents across the nodes
 
    wat_flxsa=wat_flxs;
	if(immobile) 		
		wcmo=ratio.*WC;
		wcio=WC-wcmo;
    else          
		wcmo=WC; 
		wcio=zeros(1,ncomp);
    end

    if wat_flxs(1) <=0      
		if immobile   
            wcma(1)=0.999;
		else
			wcma(1)=1.0;
        end
	else 
		if immobile   
			wcma(1)= ratio(1).*WC(1);  
		else
			wcma(1)=WC(1);
        end
    end
	if immobile   
        wcma(1,2:ncomp)=(wcmo(1,2:ncomp)+wcmo(1,1:ncomp-1))/2;
        wcma(ncomp+1)=wcma(ncomp);
    else  
		wcma(1,2:ncomp)=WC(1,2:ncomp);
		wcma(ncomp+1)=wcma(ncomp);
    end
		wco(1,1:ncomp)=WC(1,1:ncomp);

%		calculation of the pore velocity across the nodes
 
	pvela(1,1:ncomp+1)=wat_flxsa(1,1:ncomp+1)./wcma(1,1:ncomp+1);
	wat_flxso(1,1:ncomp)=(wat_flxsa(1,1:ncomp)+wat_flxsa(1,2:ncomp+1))/2; 
	pvelo(1,1:ncomp)=wat_flxso(1,1:ncomp)./wcmo(1,1:ncomp);       
    
    dm=lm.*abs(pvela');
   
    for j=1:nsol
        de(1:ncomp,j)=(dif(1:ncomp,j)).*ar(1:ncomp,j).*exp(br(1:ncomp,j).*(wcmo(1,1:ncomp))'); 
        diffus(1:ncomp,j)=de(1:ncomp,j)./(wcma(1,1:ncomp))'+dm(1:ncomp);
    end
    
end

%    calculation of the fluxes and water contents across the nodes

wat_flxsab(1,1:ncomp+1) = wat_flxsa(1,1:ncomp+1);  
wat_flxsa(1,1:ncomp+1)=   wat_flxs(1,1:ncomp+1);  
wcmob(1,1:ncomp) = wcmo(1,1:ncomp);
wciob(1,1:ncomp) = wcio(1,1:ncomp);

if immobile    
	wcmo(1,1:ncomp)=ratio(1,1:ncomp).*WC(1,1:ncomp);
	wcio(1,1:ncomp)=WC(1,1:ncomp)-wcmo(1,1:ncomp);
else
	wcmo(1,1:ncomp)= WC(1,1:ncomp);
	wcio=zeros(1,ncomp);
end
	wcmab(1,1:ncomp+1) = wcma(1,1:ncomp+1);

if wat_flxsa(1)<=0   
	if immobile 
		wcma(1)=0.999;
	else
		wcma(1)=1.0;
    end
else 
	if immobile  
		wcma(1)=ratio(1)*WC(1);  
	else
		wcma(1)=WC(1);
    end
end

if immobile 
	wcma(1,2:ncomp)=(wcmo(1,2:ncomp)+wcmo(1,1:ncomp-1))/2;
	wcma(ncomp+1)=wcma(ncomp);
else
	wcma(1,2:ncomp)=WC(1,2:ncomp);
	wcma(ncomp+1)=wcma(ncomp);
end
 
    wcob(1,1:ncomp) = wco(1,1:ncomp);
    wco(1,1:ncomp)=WC(1,1:ncomp);

%     calculation of the pore velocity across the nodes 
 
pvelab(1,1:ncomp+1) = pvela(1,1:ncomp+1);
pvela(1,1:ncomp+1)=wat_flxsa(1,1:ncomp+1)./wcma(1,1:ncomp+1); 
pvelob(1,1:ncomp) = pvelo(1,1:ncomp);
wat_flxso(1,1:ncomp)=(wat_flxsa(1,1:ncomp)+wat_flxsa(1,2:ncomp+1))/2;  
pvelo(1,1:ncomp)=wat_flxso(1,1:ncomp)./wcmo(1,1:ncomp);

%    calculation of the fluxes and pore velocity across the nodes and   half ime

wat_flxsah(1,1:ncomp+1)=(wat_flxsab(1,1:ncomp+1)+wat_flxsa(1,1:ncomp+1))/2;
wcmah(1,1:ncomp+1)=(wcmab(1,1:ncomp+1)+wcma(1,1:ncomp+1))/2;
pvelah(1,1:ncomp+1)=(pvelab(1,1:ncomp+1)+pvela(1,1:ncomp+1))/2;
pveloh(1,1:ncomp)=(pvelob(1,1:ncomp)+pvelo(1,1:ncomp))/2;

%     set the upper boundary condition
if wat_flxsa(1)>0 
	wat_flxsa(1)=0.0;
	wat_flxsab(1)=0.0;
	wat_flxsah(1)=0.0;
end
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%calculation of the apparent diffsivity and choice of the solute upper BC type

dm=lm.*abs(pvela');
%%%%%%%%%%%%%%%%%%%%
% if lm<=1.0
%    dm=max(0,(lm-0.139*dx)); % correction for numerical dispersion
% end
%dm(1:ncomp) = max(0,(lm(1:ncomp)-0.139*dx).*abs(pvela(1:ncomp)'));
%%%%%%%%%%%%%%%%%%%%
diffusb=diffus;    % apparent diffusion coefficient at t-dt 
for j=1:nsol
    de(1:ncomp,j)= (dif(1:ncomp,j)).*ar(1:ncomp,j).*exp(br(1:ncomp,j).*(wcmo(1,1:ncomp))'); % effective diffusion coefficient at time t
    diffus(1:ncomp,j)= de(1:ncomp,j)./(wcma(1,1:ncomp))'+dm(1:ncomp); 
    diffush(1:ncomp,j)= (diffusb(1:ncomp,j)+diffus(1:ncomp,j))/2;            %%%%%%%%%% apparent diffusion coefficient at half time

    [conc,solute_applic,solboco_top_type,solboco_input_type]= In_solute_boundary(t,nsol);
    if solboco_top_type == 2
    diffush(1,j)=0;   %%% set to 0 for flux BC or delete this line for concentration BC
    end
end

if  joachim == 1
    initsol = 2;
end

%% FUNCTION OUTPUT

% soil cilinder (rooted/unrooted) state params
soilStateCilinderParams.diffus = diffus;
soilStateCilinderParams.diffush = diffush;
soilStateCilinderParams.initsol = initsol;
soilStateCilinderParams.pvela = pvela;
soilStateCilinderParams.pvelo = pvelo;
soilStateCilinderParams.pvelah = pvelah;
soilStateCilinderParams.pveloh = pveloh;
soilStateCilinderParams.wat_flxsa = wat_flxsa;
soilStateCilinderParams.wat_flxsah = wat_flxsah;
soilStateCilinderParams.wcio = wcio;
soilStateCilinderParams.wciob = wciob;
soilStateCilinderParams.wcma = wcma;
soilStateCilinderParams.wcmah = wcmah;
soilStateCilinderParams.wcmob = wcmob;
soilStateCilinderParams.wcmo = wcmo;
soilStateCilinderParams.wco = wco;
soilStateCilinderParams.wcob = wcob;


