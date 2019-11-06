function [soilInnerStateParams,soilOuterStateParams] = In_initial_pressure_head(simulationSettings,...
                                                                                soilConsParams,...
                                                                                soilInnerStateParams,...
                                                                                soilOuterStateParams)
%% DOCUMENTATION                                                       
%Initial pressure head
% if i=
%		1) Initial pres. head profile calculated from top and bottom pres. head values: i=1
%		2) Initial pres. head profile calculated from some pres. head values given for some depths: i=2
%		3) Uniform initial pres. head profile: i=3
%		4) Initial pres. head profile calculated from initial water content given for some depths
%       5) Initial pres. head profile calculated from equilibrum with groundwater table
% IN:
%	compartiments_number,soil_parameters, dx_inter,BOUNDARY_CONDITIONS_MATRIX,initial_gwl 
%   dx_inter
% OUT:
%	ph: initial pressure head (L)
% CALLS:
%	none

%% FUNCTION INPUT
%simulation settings
compartiments_number = simulationSettings.nComp; % number of soil layers

%soil constant parameters
initial_gwl = soilConsParams.initial_gwl;

%soil constant parameters
ALFA = soilConsParams.alfa; % inverse air entry value [cm^-1]
N = soilConsParams.N; % van Genuchten en Mualem equation parameter
WCR = soilConsParams.wcr; % residual volumetric soil water content 
WCS = soilConsParams.wcs; % saturated volumetric soil water content




%% FUNCTION MAIN BODY
i=1;%Fill it

if i==1
   phbot = 0;	
   phtop = -85;	
   ph=phtop:(phbot-phtop)/(compartiments_number-1):phbot;
  
elseif i==2
   ph=[-369.856010443248;-368.817795024000;-367.800283699172;-366.803254734959;...
       -365.826478803138;-364.869719981445;-363.932736665022;-363.015282396446;...
       -362.117106621273;-361.237955375431;-360.377571910274;-359.535697260658;...
       -358.712070760914;-357.906430513234;-357.118513812576;-356.348057531884;...
       -355.594798471090;-354.858473673085;-354.138820709577;-353.435577939528;...
       -352.748484742620;-352.077281730020;-351.421710934508;-350.781515981874;...
       -350.156442245334;-349.546236984557;-348.950649470783;-348.369431099382;...
       -347.802335491075;-347.249118582984;-346.709538710518;-346.183356681083;...
       -345.670335840468;-345.170242132725;-344.682844154282;-344.207913202946;...
       -343.745223322442;-343.294551343030;-342.855676918733;-342.428382561644;...
       -342.012453673750;-341.607678576662;-341.213848539622;-340.830757806103;...
       -340.458203619314;-340.095986246875;-339.743909004918;-339.401778281816;...
       -339.069403561774;-338.746597448440;-338.433175688707;-338.128957196858;...
       -337.833764079180;-337.547421659162;-337.269758503388;-337.000606448204;...
       -336.739800627253;-336.487179499927;-336.242584880816;-336.005861970180;...
       -335.776859385507;-335.555429194172;-335.341426947227;-335.134711714353;...
       -334.935146119961;-334.742596380472;-334.556932342769;-334.378027523814;...
       -334.205759151419;-334.040008206178;-333.880659464519;-333.727601542878;...
       -333.580726942957;-333.439932098056;-333.305117420444;-333.176187349739;...
       -333.053050402277;-332.935619221426;-332.823810628820;-332.717545676487;...
       -332.616749699812;-332.521352371337;-332.431287755333;-332.346494363128;...
       -332.266915209156;-332.192497867688;-332.123194530225;-332.058962063512;...
       -331.999762068161;-331.945560937838;-331.896329919008;-331.852045171208;...
       -331.812687827829;-331.778244057399;-331.748705125345;-331.724067456235;...
       -331.704332696485;-331.689507777536;-331.679604979507;-331.674641995320];
   ph = ph';
   
elseif i==3
   
    ph(1:compartiments_number)=-200;
   
elseif i==4
   
   zdata=[1 85];		%depths (L) for which moisture contents are available
   wcdata=[0.3784 0.3784];    %Moisture content for each depth
   z= 1:compartiments_number;
   WC=spline(zdata,wcdata,z);%MAMADOU
   M=1-1./N;
   Se=(WC'-WCR)./(WCS-WCR);
   ph=-(Se.^(-1./M)-1).^(1./N)./ALFA;ph=ph'; %origine
   
else  %if Equilibrum with a ground water table)
    
    gwl=initial_gwl;
    ph=-abs(gwl)+abs([1:ncs].*dx_inter(1:compartiments_number));
    
end

%% FUNCTION OUTPUT
%soilInnerStateParams
soilInnerStateParams.ph = ph; % initial pressure head [cm] profile in the soil

%soilOuterStateParams
soilOuterStateParams.ph = ph; % initial pressure head [cm] profile in the soil
