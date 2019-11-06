function ph = calc_ph_from_wc(WC, soil_parameters)

WCR=soil_parameters(:,1);
WCS=soil_parameters(:,2);
ALFA=soil_parameters(:,3);
N=soil_parameters(:,4);
M=1-1./N;
Se=(WC'-WCR)./(WCS-WCR);
ph=-(Se.^(-1./M)-1).^(1./N)./ALFA;ph=ph'; %origine