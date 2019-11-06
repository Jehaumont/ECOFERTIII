function  soilCommonStateParams = calc_FC(cropStateParams,soilCommonStateParams,soilConsParams)
%% DOCUMENTATION

%% FUNCTION INPUT
% crop state parameters
drz = cropStateParams.TotRootDepth(end); % total rooting depth

% soil constant parameters
wcr    = soilConsParams.wcr; % residual volumetric water content
wcs    = soilConsParams.wcs; % saturation volumetric water content
alfa   = soilConsParams.alfa; % inverse air entry value [cm^-1]
n      = soilConsParams.N; % van Genuchten Mualem model parameter
ks     = soilConsParams.ks; % saturated hydraulic conductivity
lambda = soilConsParams.lambda; % van Genuchten Mualem model parameter
alfa_r = soilConsParams.alfa_r; % hysteris parameter of the moisture retention curve 

m      = 1-1./n; % van Genuchten Mualem parameter
ph = -100; %Field capacity is pF = 2;


%% FUNCTTION MAIN BODY
h=abs(ph.');
h(ph.'>0)=0;

WCFCi=(wcr+(wcs-wcr)./(1+(alfa.*h).^n).^m)';

WCFC = sum(WCFCi(1:drz));

%% FUNCTION OUTPUT
soilCommonStateParams.WCFC = WCFC; % volumetric water content at field capacity