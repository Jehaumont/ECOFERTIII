function [soilInnerStateParams, soilOuterStateParams] = change_state_var(fract1_new, fract2_new,...
            soilConsParams,soilInnerStateParams,soilOuterStateParams)

%% FUNCTION INPUT

soil_parameters = [soilConsParams.wcr,soilConsParams.wcs,soilConsParams.alfa,...
    soilConsParams.N,soilConsParams.ks,soilConsParams.lambda,soilConsParams.alfa_r];

wcmo1 = soilOuterStateParams.wcmo;
wcmo2 = soilInnerStateParams.wcmo;


%% FUNCTION MAIN BODY
if fract1_new==1
    ph1 = calc_ph_from_wc(wcmo1, soil_parameters);
    ph2 =[];
    
    
elseif fract2_new==1
    ph2 = calc_ph_from_wc(wcmo2, soil_parameters);
    ph1=[];
    
else
    ph1 = calc_ph_from_wc(wcmo1, soil_parameters);
    ph2 = calc_ph_from_wc(wcmo2, soil_parameters);

end

%% FUNCTION OUTPUT
soilInnerStateParams.ph = ph2;
soilOuterStateParams.ph = ph1;




