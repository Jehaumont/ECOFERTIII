function  	CH=diff_moist_capac(ph,soil_parameters)

%Calculation of the differential moisture capacity 
%
%IN:
% 	ph = soil water pressure head (cm)
%OUT:
% 	CH = differential moisture capacity (1/min)
%CALLS:
%   none
%CALLED BY:
%   state_var.m
%--------------------------
% M. Vanclooster 2/2/2000

WCR=soil_parameters(:,1);
WCS=soil_parameters(:,2);
ALFA=soil_parameters(:,3);
N=soil_parameters(:,4);

M=1-1./N;
H=abs(ph');
H(find(ph.'>=0))=1.0000e-10;					% if ponding, CH=0

% Van Genuchten model (1980), SSSAJ 44:892-898
CH=(((WCS-WCR).*M.*N.*(ALFA.*H).^N)./(H.*(1+(ALFA.*H).^N).^(M+1)))';
CH(find(isnan(CH)))=1.0000e-10;
CH(find(CH<1.0000e-10))=1.0000e-10;		%important: this line warrants that it converges
%if CH ->0, thomas block doesn't work well.