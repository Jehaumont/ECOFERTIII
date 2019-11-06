function soilStateCilinderParams = solve_flow(cropStateParams,simulationSettings,...
    managementSettings,soilConsParams,soilStateCilinderParams)

%% DOCUMENTATION
%Solve the 1-D water flow equation using finite differences
%
% IN:
%   ph = the initial pressure head (cm)
%   t = time (min)
%   dt = the time increment (min)
%   Characteristics of the soil: ncomp,dx,dx_inter,soil_parameters,stock_max
%   General characteristics:maxiter,dt_min,err_tol,dt_start,hs,phsa,err_tol,simplant
%   Current boundary specifications:boco_top_type,boco_top,boco_bot_type,boco_bot,seep,ponded,
%       pond_from,iter,bctop_changed,bcbot_changed
%   Current soil values: phsurf,phbot,flxsbot,flxar,pond,flxa1
% OUT:boco_top_type
%   ph = the updated pressure head (cm)
%   wat_flxs = flow at each node (cm/min)
%   dt = the update time step (min)
%   iter = number of solve_flow iterations
% CALL:
%   state_var,conduct_in,thomas_block2,calc_fluxes,check_balance,calc_dt
% CALLED BY:
%   wavemat101.m
%
%----------------------------------
% M. Vanclooster, 13/1/2000
% modified by M. Javaux, 14/05/00, updating 1:17-11-00
% modified by M.Sall, 25/11/09

%% FUNCTION INPUT
% crop state params
drz = cropStateParams.TotRootDepth(end);

% simulation settings
t = simulationSettings.t;
ncomp = simulationSettings.nComp;
dt_min = simulationSettings.dt_min;
dx = simulationSettings.dx;
dx_inter = simulationSettings.dxInter;
err_tol = simulationSettings.errTol;
maxiter = simulationSettings.maxIter;
simplant = simulationSettings.simPlantFlag;
units = simulationSettings.units;

%management settings

harvest_date = sort([managementSettings.hDateCauli;managementSettings.hDateLeek]);
plant_date = sort([managementSettings.pDateCauli;managementSettings.pDateLeek]);

% soil constant parameters
hs = soilConsParams.hs;
phsa = soilConsParams.phsa;
soil_parameters = [soilConsParams.wcr,soilConsParams.wcs,soilConsParams.alfa,...
    soilConsParams.N,soilConsParams.ks,soilConsParams.lambda,soilConsParams.alfa_r];

% soil state cilinder params
boco_bot = soilStateCilinderParams.boco_bot;
boco_bot_type = soilStateCilinderParams.boco_bot_type;
boco_top = soilStateCilinderParams.boco_top;
boco_top_type = soilStateCilinderParams.boco_top_type;
dt = soilStateCilinderParams.dt;
epa = soilStateCilinderParams.epa;
flxa1 = soilStateCilinderParams.flxa1;
flxsbot = soilStateCilinderParams.flxsbot;
flxar = soilStateCilinderParams.flxar;
iter = soilStateCilinderParams.iter;
ph = soilStateCilinderParams.ph;
phbot = soilStateCilinderParams.phbot;
phsurf = soilStateCilinderParams.phsurf;
pond = soilStateCilinderParams.pond;
pond_from = soilStateCilinderParams.pond_from;
ponded = soilStateCilinderParams.ponded;
stock_max = soilStateCilinderParams.stock_max;

%initialization
%put the "changement variables" equal to zero
bcbot_changed=0;
bctop_changed=0;
DX(1:ncomp)=dx;
iter_dx=0;
itertop=0;
top_OK=0;
bot_OK=0;
no_conv=0;

%% FUNCTION MAIN BODY
%initial calculations
[WC,kh,CH,rtex,EPRA]=state_var(cropStateParams,managementSettings,simulationSettings,...
    soilConsParams,soilStateCilinderParams,0);
phB = ph;
WCB = WC;
kh_in=conduct_in(ph,boco_top_type,phsurf,phsa,flxar,soil_parameters);


%calculates current boundary conditions
[seep,boco_top_type,boco_bot_type,dt,bctop_changed,...
    bcbot_changed,flxar,phbot,phsurf,flxsbot]=calc_boco(dx,kh_in,ph,dt,rtex,t,...
    ncomp,pond,boco_bot,boco_top_type,boco_bot_type,...
    flxa1,stock_max,bctop_changed,bcbot_changed,soil_parameters,flxar,...
    phsa,phsurf,phbot,flxsbot,ponded,pond_from);

%Resolves Thomas block by Newton-Raphson
temp = dt;
if temp < dt_min
    dt = dt_min;
end

while (dt >= dt_min) && (top_OK==0) && (itertop<=3)
    dt = temp;
    if itertop==2
        [phsurf,boco_top_type]=fix_uboco(phsa,pond,flxar);
        disp('fix')
    end
    itertop=itertop+1;
    ph = phB;
    soilStateCilinderParams.ph = ph;
    [WC,kh,CH,rtex,EPRA]=state_var(cropStateParams,managementSettings,simulationSettings,...
        soilConsParams,soilStateCilinderParams,0);
    kh_in=conduct_in(ph,boco_top_type,phsurf,phsa,flxar,soil_parameters);
    
    % in case of free drainage,re(calculate) flxsbot
    if boco_bot_type==4
        flxsbot=-kh_in(ncomp+1);
    end
    balance_error = 1;
    
    while balance_error  == 1
        iter = 0;
        
        while balance_error==1 && (iter < maxiter)
            iter = iter + 1;
            % call thomasblock to inverse tridiag. matrix
            ph = thomas_block2(ph,WC,WCB,kh,kh_in,CH,rtex,dt,DX,dx_inter,...
                phsurf,phbot,flxsbot,flxar,boco_top_type,boco_bot_type,...
                ncomp);
            
            %recalculate state variables: kh,WC,...
            soilStateCilinderParams.ph = ph;
            [WC,kh,CH,rtex,EPRA]=state_var(cropStateParams,managementSettings,simulationSettings,...
                soilConsParams,soilStateCilinderParams,0);
            %Calculation of the hydraulic conductivity in between the nodes
            kh_in=conduct_in(ph,boco_top_type,phsurf,phsa,flxar,soil_parameters);
            %Calculation of the soil water fluxes across the soil nodes
            [wat_flxs,iter] = calc_fluxes (ph,kh_in,dt,phsurf,pond,flxsbot,phbot,flxar,...
                boco_top_type,boco_bot_type,maxiter,dx_inter,ncomp,iter);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ncomp~=ncomp
                for k=ncomp:ncomp
                    WC(k)=moist_ret(0,soil_parameters(k,:),dt,0);
                    %          wat_flxs(k+1)=wat_flxs(k)+((WC(k)-WCB(k))/dt+rtex(k))*dx;
                    wat_flxs(k+1)= ((WC(k)-WCB(k))./dt+rtex(k))*dx +wat_flxs(k);
                    ph(k)=ph(k-1)+(wat_flxs(k)/kh_in(k)+1.)*dx_inter(k);
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %check of the balance
            balance_error = check_balance(wat_flxs,WC,WCB,rtex,dt,dx,ncomp,err_tol);
            % in case of free drainage,re(calculate) flxsbot
            if boco_bot_type==4
                flxsbot=-kh_in(ncomp+1);
                
            end
        end
        
        if iter == maxiter && balance_error
            dt=dt/2;
            if dt < dt_min
                no_conv=1;
                break;
            else
                disp('time step diminution');
                ph = phB ;
                soilStateCilinderParams.ph = ph;
                [WC,kh,CH,rtex,EPRA] = state_var(cropStateParams,managementSettings,simulationSettings,...
                                        soilConsParams,soilStateCilinderParams,0);
                kh_in = conduct_in(ph,boco_top_type,phsurf,phsa,flxar,soil_parameters);
                [seep, boco_top_type, boco_bot_type, dt, bctop_changed,...
                    bcbot_changed, flxar, phbot, phsurf, flxsbot] = ...
                    calc_boco(dx,kh_in,ph,dt,rtex,t,ncomp,...
                    pond,boco_bot,boco_top_type,boco_bot_type,...
                    flxa1,stock_max,bctop_changed,bcbot_changed,...
                    soil_parameters,flxar, phsa,phsurf,phbot,flxsbot,ponded,pond_from);
                
            end
        end
    end
    [top_OK, bot_OK, bctop_changed, bcbot_changed, phsurf, phbot, flxsbot,boco_top_type,...
        boco_bot_type,boco_top]=check_bc(wat_flxs,ph,kh_in,ncomp,dx_inter,pond,phsa,...
        flxar,boco_top_type,boco_bot_type,seep,phsurf,phbot,flxsbot,bctop_changed,...
        bcbot_changed,boco_top);
end

%ponding
pond = min([hs  max([0 (wat_flxs(1)-flxar)*dt])]);

%in case of hysteresis
if soil_parameters(1,7)~=1
    WC = moist_ret(ph, soil_parameters, dt, 1);
end

%Runoff generation
if pond <hs && pond >0
    disp(sprintf('Ponding:%4.5f cm',pond));
end

if (pond==hs & pond~=0)|(hs==0 & (wat_flxs(1)-flxar)*dt>0)
    runoff=((wat_flxs(1)-flxar)*dt-hs)/dt;
    disp(['Runoff generation: ',num2str(runoff),' cm/minutes (=',num2str(abs(runoff/flxar)),'% of the prescribed flux)']);
else
    runoff=0;
end

%% FUNCTION OUTPUT

% soil state cilinder params
soilStateCilinderParams.bctop_changed = bctop_changed;
soilStateCilinderParams.bcbot_changed = bcbot_changed;
soilStateCilinderParams.boco_bot_type = boco_bot_type;
soilStateCilinderParams.boco_top = boco_top;
soilStateCilinderParams.boco_top_type = boco_top_type;
soilStateCilinderParams.dt = dt;
soilStateCilinderParams.EPRA = EPRA;
soilStateCilinderParams.iter = iter;
soilStateCilinderParams.flxar = flxar;
soilStateCilinderParams.no_conv = no_conv;
soilStateCilinderParams.ph = ph;
soilStateCilinderParams.phsurf = phsurf;
soilStateCilinderParams.pond = pond;
soilStateCilinderParams.rtex = rtex;
soilStateCilinderParams.runoff = runoff;
soilStateCilinderParams.wat_flxs = wat_flxs;
soilStateCilinderParams.WC = WC;

