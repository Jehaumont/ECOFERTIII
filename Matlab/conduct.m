function  kh = conduct(ph,soil_parameters)

%Calculation of the hydraulic conductivity
%IN:
% 	ph = soil water pressure head (m)
%OUT:
% 	kh = soil hydraulic conductivity  (m/s)
%CALL:
%   none
%CALLED BY:
%   conduct_in.m
%----------------------
% M. Vanclooster 2/2/2000

WCR=soil_parameters(:,1);
WCS=soil_parameters(:,2);
ALFA=soil_parameters(:,3);
N=soil_parameters(:,4);
KS=soil_parameters(:,5);
LAMBDA=soil_parameters(:,6);
M=1-1./N;
H=abs(ph');
H(find(ph.'>0))=0.0;
% Van Genuchten model (1980), SSSAJ 44:892-898
SE=(1+(ALFA.*H).^N).^(-M);
kh=(KS.*SE.^LAMBDA.*(1-(1-SE.^(1./M)).^M).^2)';
kh(find(isnan(kh)))=1.0000e-10;
		%kh(find(kh<1.0000e-10))=1.0000e-10;%changed from 
kh(find(kh==0))=1.0000e-10; 
