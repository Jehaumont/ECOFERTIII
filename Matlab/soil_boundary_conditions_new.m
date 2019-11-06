function soilStateCilinderParams=soil_boundary_conditions_new(climateState,simulationSettings,soilConsParams,soilStateCilinderParams)
%% DOCUMENTATION
% Calculate/find the current boundary conditions data
%
% IN :
%   General characteristics:t,ph,dt,pot_surface_flux,simplant
%   Characteristics of the bc in the last iteration: boco_top_type,boco_top,boco_bot_type,boco_bot,pond,runoff,phsurf
%   Changes during last iteration: var_changed (=bctop_changed,bcbot_changed,dt_changed_bc,first_time_bc)
% OUT :
%   dt,hs,phsa,pot_surface_flux,flxa1,boco_top_type,boco_top,boco_bot_type,
%   boco_bot,pond,bctop_changed,bcbot_changed,first_time_bc,dt_changed
% CALLS:
%   In_Boundary_conditions.m,In_ETsplit.m
% CALLED BY:
%   wavemat105.m
%------------------------------------------------
% Javaux M., Lambot S. & Vanclooster M. (2008)
%modified by M.Sall 25/11/09

%% FUNCTION INPUT

% simulation settings
t = simulationSettings.t;
dx_inter = simulationSettings.dxInter;
ncomp = simulationSettings.nComp;

% soil common state params
arel = soilConsParams.arel;
brel = soilConsParams.brel;

% soil state cilinder (root/unrooted) parameters
bctop_changed = soilStateCilinderParams.bctop_changed;
boco_top = soilStateCilinderParams.boco_top;
boco_top_type = soilStateCilinderParams.boco_top_type;
BOUNDARY_CONDITIONS_MATRIX = soilStateCilinderParams.BOUNDARY_CONDITIONS_MATRIX;
dt = soilStateCilinderParams.dt;
dt_changed_bc  = soilStateCilinderParams.dt_changed_bc;
first_time_bc = soilStateCilinderParams.first_time_bc;
next_dt_new_bc = soilStateCilinderParams.next_dt_new_bc;
ph = soilStateCilinderParams.ph;
phsurf = soilStateCilinderParams.phsurf;
pond = soilStateCilinderParams.pond;
runoff  = soilStateCilinderParams.runoff;

%initialization
if first_time_bc
    dt_changed_bc=0;
    first_time_bc=0;
    next_dt_new_bc=0;
end

info_pond=0;

% Searches for the appropriate boundary condition for the current time t 

i=max(find(BOUNDARY_CONDITIONS_MATRIX(:,1)<=t));

if dt_changed_bc==1 % In case the dt have been changed to fit the limits of the BC (see calcdt)
    dt_changed_bc=0;    
    if next_dt_new_bc==1
        next_dt_new_bc=0;
        [~,soilStateCilinderParams] = In_Boundary_conditions(1,climateState,soilConsParams,soilStateCilinderParams); %TO DO find a way to check original input of BC
        if(soilStateCilinderParams.IPBCM(i,2)==2 && boco_top_type==1 && boco_top>0)
        %if (BOUNDARY_CONDITIONS_MATRIX(i,2)==2 & boco_top_type==1 & boco_top>0)
        % New user BC is now flux (2) and  different from previous (1)
        % and there is water which has not been infiltrated (boco_top>0);
        % In this case, ponding will appear. 
        % In all other cases, no recalculations are needed
            imposed_flux=BOUNDARY_CONDITIONS_MATRIX(i,3);
            info_pond=1;
        %next BC will be ponding
            BOUNDARY_CONDITIONS_MATRIX(i,2)=1;
            BOUNDARY_CONDITIONS_MATRIX(i,3)= boco_top; % previous reserve is still considered     
        end
    end
    next_dt_new_bc=1;
end    

[~,soilStateCilinderParams] = In_Boundary_conditions(1,climateState,soilConsParams,soilStateCilinderParams); %TO DO find a way to check original input of BC
%if IPBCM(i,2) ~= BOUNDARY_CONDITIONS_MATRIX(i,2)
if bctop_changed==1
% Current BC are different of user BC
% The calculated flux by solving the flow equation (calc_boco, check_bc) 
% is incompatible with imposed BC
if soilStateCilinderParams.IPBCM(i,2) ~=boco_top_type
    if soilStateCilinderParams.IPBCM(i,2)==2
    %if BOUNDARY_CONDITIONS_MATRIX(i,2)==2
        % if user BC is flux (current is ph), then define the imposed flux
        imposed_flux=BOUNDARY_CONDITIONS_MATRIX(i,3);
        %imposed_flux=IPBCM(i,3);
        info_pond=1;
    end
    %definition of the future BC=current BC
    BOUNDARY_CONDITIONS_MATRIX(i,2)=boco_top_type;
    %boco_top_type = BOUNDARY_CONDITIONS_MATRIX(i,2);
    if boco_top_type==1									
    % new TBC is pression/succion
        BOUNDARY_CONDITIONS_MATRIX(i,3)=phsurf;
    else														
         % new TBC is flux
      	%BOUNDARY_CONDITIONS_MATRIX(i,3)=pot_surface_flux;  
    end
end
end

%Next Boundary condition for time t
boco_top_type=BOUNDARY_CONDITIONS_MATRIX(i,2) ;
boco_top=BOUNDARY_CONDITIONS_MATRIX(i,3);
boco_bot_type=BOUNDARY_CONDITIONS_MATRIX(i,4);
boco_bot=BOUNDARY_CONDITIONS_MATRIX(i,5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% In case of the presence of a ground water table
A= BOUNDARY_CONDITIONS_MATRIX(i,4);
B= BOUNDARY_CONDITIONS_MATRIX(i,5);
if A== 5 || A== 6 || A== 7
 [~,boco_bot,~,BOUNDARY_CONDITIONS_MATRIX]=...
    find_gwl(A,B,ph,dx_inter,ncomp,BOUNDARY_CONDITIONS_MATRIX,phsurf,arel,brel);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Upper boundary condition
if boco_top_type==1
   %succion/pression
    [~,soilStateCilinderParams] = In_Boundary_conditions(1,climateState,soilConsParams,soilStateCilinderParams); %TO DO find a way to check original input of BC
   if soilStateCilinderParams.IPBCM(i,2)==2
   	imposed_flux = soilStateCilinderParams.IPBCM(i,3);
   else
       imposed_flux =0;
   end
   flxa1=imposed_flux;%esa+imposed_flux;
   phsurf=boco_top;
elseif boco_top_type==2										
   % flux, pluie
   flxa1=boco_top;%esa+boco_top;
end

if info_pond==1
   %ponding appears
   flxa1=imposed_flux;%esa+imposed_flux;
end

%potential flux on the top (=flxar)
pot_surface_flux=flxa1-pond/dt;
%that is the water which is imposed by the user and the water in the ponding

%if runoff, pot surf flux is reducted
if runoff~=0
   pot_surface_flux=pot_surface_flux-runoff;
end

%% FUNCTION OUTPUT

soilStateCilinderParams.boco_bot = boco_bot;
soilStateCilinderParams.boco_bot_type = boco_bot_type;
soilStateCilinderParams.BOUNDARY_CONDITIONS_MATRIX = BOUNDARY_CONDITIONS_MATRIX;
soilStateCilinderParams.flxar = pot_surface_flux;
soilStateCilinderParams.flxa1 = flxa1;
soilStateCilinderParams.boco_top_type = boco_top_type;
soilStateCilinderParams.boco_top = boco_top;
soilStateCilinderParams.pond = pond;
soilStateCilinderParams.bctop_changed  = bctop_changed;
soilStateCilinderParams.first_time_bc = first_time_bc;
soilStateCilinderParams.dt_changed_bc = dt_changed_bc;
soilStateCilinderParams.next_dt_new_bc = next_dt_new_bc;
soilStateCilinderParams.phsurf = phsurf;

