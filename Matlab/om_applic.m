function om_appl = om_applic(managementSettings,soilStateCilinderParams)

%% FUNCTION INPUT

hDateCauli = managementSettings.hDateCauli;
om_appl = soilStateCilinderParams.om_appl;

%% FUNCTION MAIN BODY

om_appl = [hDateCauli,repmat(om_appl,length(hDateCauli),1)];
om_appl = [zeros(1,size(om_appl,2));om_appl];


