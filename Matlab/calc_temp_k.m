function aka = calc_temp_k(geomq, geomom,geomos, alamo,alama,WC,geoma);

%Calculate k coefficients for quartz, organic matter, other solids and air

%IN:
% geomq     = g-factors quartz(-)
% geomom    =  g-factors organic matter (-)
% geomos    = g-factors other solids(-)
% alamo     = thermal conductivities ofthe fractions(J L-1 T-1 C-1)
% alama     = thermal conductivity of the medium(J L-1 T-1 C-1)
% WC        = water content of the current time step(L3 L-3)
% geoma     = g-factors air(-)

%OUT:
%aka : constants as defined in equation 4-4

%CALLED BY:
%calc_thermal_conduc

for m=1:3;
    if m==1
        g = geomq;
    elseif m==2
        g=geomom;
    elseif m==3
        g=geomos;
    end
    sums = 1./(1+((alamo(:,m)./alama) -1).*g(:,1))+...
            1./(1+((alamo(:,m)./alama) -1).*g(:,2))+...
            1./(1+((alamo(:,m)./alama) -1).*g(:,3));        
    aka(:,m)=1/3*sums;
end

%Calculate lambda_v
dummy = WC> 0.2;
alamv = dummy*alamo(5) +(1-dummy).*(WC*alamo(5)/0.2);
     %decrease from 0.176 mCal cm-1 sec-1 °C to 0 when WC drops from 0.2 to 0
     %in wave fortran different: WC(i)*0.2/alamo(5)
%Fortran: alamv = dummy*alamo(5) +(1-dummy).*(WC*0.2/alamo(5));

%Calculate k-coefficient for air  
    sums = 1./(1+(((alamo(:,4)+alamv)./alama) -1).*geoma(:,1))+...
            1./(1+(((alamo(:,4)+alamv)./alama) -1).*geoma(:,2))+...
            1./(1+(((alamo(:,4)+alamv)./alama) -1).*geoma(:,3));        
    aka(:,4)=1/3*sums;
    
  