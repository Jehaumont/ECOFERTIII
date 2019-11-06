function geoma = calc_temp_g_air(WC, wcs)

%Calculate g-factors of air

%IN:
% WC: water content of the current time step (L3 L-3)
% wcs: saturated water content (L3 L-3)
%OUT:
% geoma : g-factor air (-) 
% CALLS: none 
% CALLED BY:calc_thermal_conduc


dummy = WC> 0.2;
geoma(:,1)=dummy.*((WC-wcs)*(0.105-0.333)./(0.2-wcs) +0.333) +...
    (1-dummy).*((WC-0.2)*(0.015-0.105)/(0-0.2) +0.105);
        %decrease from 0.333 to 0.105 when WC drops from saturated to 0.2
        %FORTRAN % geoma(1)= (WC(i) -0.0845)/1.10; 
        %decrease from 0.105 to 0.015 when WC drops from 0.2 to 0
        %FORTRAN: geoma(1) = (WC(i) - 0.0333)/2.22;
        
%Fortran:
% geoma(:,1) = dummy.*(WC-0.0845)/1.10 +...
%     (1-dummy).*((WC -0.0333)/2.22) ;        
 geoma(:,2) = geoma(:,1);
 geoma(:,3) = 1- geoma(:,1) - geoma(:,2);