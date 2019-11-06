function [allam, beta] = calc_thermal_conduc(units,ncomp,soil_parameters, solute_param1,WC,wcp);

%Calculates thermal conductivity

%In :
% units: specify which units are used in in_general_data
% ncomp: number of compartments
% soil_parameters(:,2) :saturated water conten
% solute_param1(:,1) : bulk density
% WC: water content of the current time step
% wcp: water content of the previous time step

%OUT:
% allam: thermal conductivity between the nodes, halfway the time step 
%in (J L-1 T-1 C-1)
%beta : heat capacity of the soil (J/C/L^3)

% CALLS:
% In_temperature_param
% calc_temp_k
% calc_temp_g_air
% CALLED BY:
% solve_temperature_JV


%From other modules:
wcs = soil_parameters(:,2); %saturated water content
bd = solute_param1(:,1); %bulk density
bd(ncomp+1) = bd(ncomp);
WC = WC'; %Water content

%Define extra water contents for calculation allam (ncomp +1) 
%Joachim V. (Not included in fortran)
WC(ncomp +1)= WC(ncomp);
wcp(ncomp +1) = wcp(ncomp);
wcs(ncomp+1) = wcs(ncomp);


    %%calculate water content halfway the node, halfway the time step
    i = [2:ncomp+1];
    tw(1) = (WC(1) + wcp(1))/2;   
    tw(i) = (WC(i) + wcp(i) + WC(i-1) + wcp(i-1))/4;
    tw = tw';

%User defined input of temperature parameters
[part, alamo, geomq, geomos,geomom,bd_water,heatcapsolids_corr, heatcapwater_corr]= In_temperature_param(units,ncomp);    
    
%Calculate volumetric heat capacity (eq 4-2 Wave manual)
i = [1:ncomp+1];
beta(i) = heatcapsolids_corr*bd(i) + heatcapwater_corr*tw(i)*bd_water; %in J/C/L^3 
beta = beta';

%calculate lambda 0 (thermal conductivity of the medium 
%(air or water with threshold set to WC = 0.2)
% value for water if > 0.2, air if < 0.2
dummy = WC> 0.2;
alama = dummy*alamo(1,5) +(1-dummy)*alamo(1,4);
alama = dummy*alamo(1,5) +(1-dummy)*alamo(1,5);

%Calculate g-factors of air
geoma = calc_temp_g_air(WC, wcs);

%Parameters needed in eq 4-3
%k coefficients for quartz, organic matter,other solids and air
aka = calc_temp_k(geomq, geomom,geomos, alamo,alama,WC,geoma);

%Fractions
%Input as variable part for 3 species
%air and medium fractions need to be calculated from WC
%Air
part(:,4) = 1- WC - part(:,1) -part(:,2) -part(:,3);

%Medium
dummy = WC> 0.2;
fx = dummy.*tw + (1-dummy).*(1- tw - part(:,1) -part(:,2) - part(:,3));

%Fill in equation 4-3
%numerator
a12 = sum(aka(1:end,1:4).*part(1:end,1:4).*alamo(1:end,1:4),2)+alama.*fx;
%denominator
a13 = sum(aka(:,1:4).*part(:,1:4),2)+fx;

%thermal conductivity
dummy = WC> 0.2;
%c = 1.65; %Correction factor
c=1;
allam = dummy.*(a12./a13)+c*(1-dummy).*(a12./a13);

