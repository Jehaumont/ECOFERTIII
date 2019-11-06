function  wc = moist_ret(ph,soilStateCilinderParams,soilConsParams,update_historic)

%% DOCUMENTATION
%Calculation of the moisture retention
%
%In:
% 	ph = soil water pressure head (L)
%	soil_parameters
%   dt:time step (T)
%   update_historic: =1 if 1st call in the iteration and 0 if not
% OUT:
% 	wc = volumetric soil moisture content
%CALLS:
%   none
%CALLED BY:
%   state_var.m
%   calc_stock.m
%-----------------------------
%by Javaux m., Vanclooster m., hysteresis by Lambot S.

%% FUNCTION INPUT
% soil cilinder state
dt = soilStateCilinderParams.dt;

% soil constant parameters
alfa = soilConsParams.alfa; %inverse air entry value [cm^-1]
alfa_r = soilConsParams.alfa_r; %hysteris parameter of moisture retention curve [-]
ks = soilConsParams.ks; % saturated hydraylic conductivity
lambda = soilConsParams.lambda; %van Genuchten Mualem model parameter
n = soilConsParams.N; %van Genuchten Mualem modle parameter
m = 1-1./n; %van Genuchten Mualem modle parameter
tolerance = soilConsParams.tolerance; % tolerance on state difference for detection fo drying or wetting
wcr = soilConsParams.wcr; % residual volumetric soil water content
wcs = soilConsParams.wcs; % saturated volumetric soil water content

%% FUNCTION MAIN BODY
%hysteresis part
if alfa_r~=1
    % Detection of state change (drying or wetting)
    if isempty(phprec)
        state=zeros(1,length(ph));							%initial condition : drying
        update_historic=1;
    else
        state=stateprec;
        if update_historic==1
            state(find((ph-phprec)/dt>+tolerance))=1;	%wetting
            state(find((ph-phprec)/dt<-tolerance))=0;	%drying
        end
    end
    % Reversal point values
    if isempty(phrevers)
        phrevers=-ones(1,length(ph));						%origin of the main drying curve
        wcrevers=wcs;											%
    else
        phrevers(find(state~=stateprec))=phprec(find(state~=stateprec));
        wcrevers(find(state~=stateprec))=wcprec(find(state~=stateprec));
    end
    
    % Parameters transformation for hysteresis
    alfa = (alfa_r.*alfa).*(state==1)+ (alfa).*(state==0);
    
    h = phrevers;
    h(find(h>0)) = 0;
    h = abs(h);
    se = (1+(alfa.*h).^n).^(-m);
    se(find(se>0.999))=0.9999;	%approximation to avoid divide by zero
    se(find(se<0.001))=0.0001;	%approximation to avoid divide by zero
    
    wcr = (wcr).*(state==0) + ((wcrevers-wcs.*se)./(1-se)).*(state==1);
    wcs = (wcs).*(state==1) + ((wcrevers-wcr.*(1-se))./se).*(state==0);
end

% Van Genuchten model (1980), SSSAJ 44:892-898
h=abs(ph.');
h(find(ph.'>0))=0;

wc=(wcr+(wcs-wcr)./(1+(alfa.*h).^n).^m)';
wc(find(isnan(wc)))=1.0000e-10;
wc(find(wc==0))=1.0000e-10;

if (update_historic)&&(alfa_r~=1)
    stateprec=state;
    phprec=ph;
    wcprec=wc;
end
%% FUNCTION OUTPUT
% soil cilinder state parameters
wc = wc;
