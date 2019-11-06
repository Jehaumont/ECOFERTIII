function soilStateCilinderParams = store_data(climateState,simulationSettings,soilStateCilinderParams)
%% DOCUMENTATION
%Store the output of solve_flow for each time step
%IN:
% 	t,ph,wc_profile,tnode,ph_node,wc_node,print_time,
%	print_node,cum_bot_flxs,bot_flux,no_conv,soil_parameters
%OUT:
%	wc_profile,tnode,wc_node,ph_node,bot_flux
%CALLS:
% moist_ret
% if no convergence, then  output are not stored
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% FUNCTION INPUT
% climate state
climate = climateState.climateDaily;

% simulation settings
dx = simulationSettings.dx;
print_node = simulationSettings.print_node;
print_time = simulationSettings.print_time;
t = simulationSettings.t;

% soil state cilinder parameters
bot_flux = soilStateCilinderParams.bot_flux;
bot_inf  = soilStateCilinderParams.bot_inf;
cum_bot_flxs = soilStateCilinderParams.cum_bot_flxs;
cum_evap = soilStateCilinderParams.cum_evap;
cum_infiltr = soilStateCilinderParams.cum_infiltr;
cum_pot_transp = soilStateCilinderParams.cum_pot_transp;
cum_potential_surface_flux  = soilStateCilinderParams.cum_potential_surface_flux;
cum_sink_wat = soilStateCilinderParams.cum_sink_wat;
cum_sink_wat2 = soilStateCilinderParams.cum_sink_wat2;
cum_top_flxs = soilStateCilinderParams.cum_top_flxs;
cum_trans = soilStateCilinderParams.cum_trans;
dt = soilStateCilinderParams.dt;
epa = soilStateCilinderParams.epa;
esa = soilStateCilinderParams.esa;
evap = soilStateCilinderParams.evap;
flxar = soilStateCilinderParams.flxar;
no_conv = soilStateCilinderParams.no_conv;
ph = soilStateCilinderParams.ph;
ph_node = soilStateCilinderParams.ph_node;
ph_profile = soilStateCilinderParams.ph_profile;
potential_surface_flux = soilStateCilinderParams.potential_surface_flux;
potential_transp  = soilStateCilinderParams.potential_transp;
sink = soilStateCilinderParams.sink;
sink_prof = soilStateCilinderParams.sink_prof;
snode = soilStateCilinderParams.snode;
tnode = soilStateCilinderParams.tnode;
top_flux = soilStateCilinderParams.top_flux;
top_inf = soilStateCilinderParams.top_inf;
trans = soilStateCilinderParams.trans;
wat_flxs = soilStateCilinderParams.wat_flxs;
water_storage = soilStateCilinderParams.water_storage;
WC = soilStateCilinderParams.WC;
wc_node = soilStateCilinderParams.wc_node;
wc_profile = soilStateCilinderParams.wc_profile;

%% FUNCTION MAIN BODY

if no_conv ~= 1
    
    %wc_profile information
    for i=1:length(print_time)
        if (t>print_time(i)) & (wc_profile(i,:)==1)
            wc_profile(i,:) = WC;
            ph_profile(i,:) = ph;
            sink(i,:) = sink_prof;
        end
    end
    
    %Observation nodes
    tnode(end+1) = t;
    wc_node(:,end+1)= (WC(print_node)).';
    ph_node(:,end+1)= (ph(print_node)).';
    
	if isempty(snode)
        snode(:,end+1) = (cum_sink_wat2(print_node)).';
    else
        snode(:,end+1) = (cum_sink_wat2(print_node)).' - snode(:,end);
    end
    
    %Cumulative bottom flux
    bot_flux(end+1) = cum_bot_flxs;
    
    %Instanate bottom flux
    bot_inf(end+1) = wat_flxs(end);
    
    %Cumulative top flux = cumulative actual surface flux
    top_flux(end+1)= -cum_top_flxs;
    
    %INstantanate top flux = Actual surface flux
    top_inf(end+1)=  -wat_flxs(1);
    
    %Potential surface flux
    potential_surface_flux(end+1) = -flxar;
    
    %Cumulative potential surface flux
    cum_potential_surface_flux(end+1) = cum_potential_surface_flux(end) + flxar*dt;
    
    %Potential transpiration
    potential_transp(end+1) = epa;

    %Cumulative potential transpiration
    cum_pot_transp(end+1) = cum_pot_transp(end) + epa*dt;
    
    %Actual transpiration;
    trans(end+1) = sink_prof(end)*dx;
    
    %cumulative actual transpiration
    cum_trans(end+1)= sum(cum_sink_wat)*dx;
    
    %Potential evaporation
    evap(end+1)=esa;
    
    %Potential infiltration (=rinf)(not for output, but calculating cum
    %infiltra and evap)
    i = max(find(climate(:,1)<=t));
    rinf = (climate(i,3) + climate(i,4));
    
    %Cumulative actual infiltration and actual evaporation
    if wat_flxs(1) > 0 && esa ~= 0
        
		cum_evap(end+1) = cum_evap(end) + abs(wat_flxs(1)) * dt;
        cum_infiltr(end+1) = cum_infiltr(end);
        
		if rinf > 0
            cum_evap(end) = cum_evap(end) + rinf * dt;
            cum_infiltr(end) = cum_infiltr(end) + rinf * dt;
        end
    else
        
		cum_infiltr(end+1) = cum_infiltr(end) + abs(wat_flxs(1)) * dt + esa * dt;
        cum_evap(end+1) = cum_evap(end)+ esa * dt;
    end
    
else
    no_conv = 0;
end

%Water storage in the profile
water_storage(end+1) = sum(WC.*dx);

%% FUNCTION OUTPUT

% soil state cilinder parameters
soilStateCilinderParams.wc_profile = wc_profile;
soilStateCilinderParams.tnode = tnode;
soilStateCilinderParams.wc_node = wc_node;
soilStateCilinderParams.ph_node = ph_node;
soilStateCilinderParams.bot_flux = bot_flux;
soilStateCilinderParams.top_flux = top_flux;
soilStateCilinderParams.top_inf = top_inf;
soilStateCilinderParams.snode = snode;
soilStateCilinderParams.evap = evap;
soilStateCilinderParams.trans = trans;
soilStateCilinderParams.sink = sink;
soilStateCilinderParams.ph_profile = ph_profile;
soilStateCilinderParams.bot_inf  = bot_inf ;
soilStateCilinderParams.potential_surface_flux  = potential_surface_flux ;
soilStateCilinderParams.cum_potential_surface_flux  = cum_potential_surface_flux ;
soilStateCilinderParams.potential_transp  = potential_transp ;
soilStateCilinderParams.cum_evap  = cum_evap ;
soilStateCilinderParams.cum_infiltr = cum_infiltr;
soilStateCilinderParams.cum_pot_transp = cum_pot_transp;
soilStateCilinderParams.cum_trans = cum_trans;
soilStateCilinderParams.water_storage = water_storage;
