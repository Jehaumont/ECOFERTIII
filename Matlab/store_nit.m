function soilStateCilinderParams = store_nit(soilStateCilinderParams)

%% FUNCTION INPUT
cberr_corg = soilStateCilinderParams.cberr_corg;
cberr_norg = soilStateCilinderParams.cberr_norg;
cum_nit_sink = soilStateCilinderParams.cum_nit_sink;
ptden = soilStateCilinderParams.ptden;
pthyd = soilStateCilinderParams.pthyd;
ptmin = soilStateCilinderParams.ptmin;
ptnit = soilStateCilinderParams.ptnit;
ptscorg = soilStateCilinderParams.ptscorg;
ptsnorg  = soilStateCilinderParams.ptsnorg; 
ptup = soilStateCilinderParams.ptup;
ptvol  = soilStateCilinderParams.ptvol;
tflcorg = soilStateCilinderParams.tflcorg;
tflnorg = soilStateCilinderParams.tflnorg;

%% FUNCTION MAIN BODY
%Store nitrogen sink terms
cum_nit_sink(end+1,:) = [ptsnorg, ptscorg, ptup', ptmin', pthyd,ptnit,ptvol,ptden, tflnorg, tflcorg, cberr_norg, cberr_corg];

%% FUNCTION OUTPUT
soilStateCilinderParams.cum_nit_sink = cum_nit_sink;