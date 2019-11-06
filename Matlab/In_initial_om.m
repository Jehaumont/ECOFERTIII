function [soil_om]=In_initial_om(soil_parameters,ncs);

%Initial organic matter content of the profile

%IN
%soil_parameters, ncs(number of compartiments)
%OUT
%C_MANP: carbone in manure (M L-2)
%N_MANP: nitrogen in manure (M L-2)
%C_LITP: carbone in litter (M L-2)
%N_LITP: nitrogen in litter (M L-2)
%C_HUMP: carbone in manure (M L-2)
%N_HUMP: nitrogen in humus (M L-2)

%CALL:none
%CALLED BY: WAVE_MAT
%--------------------------------------------------------------------------
%M.SALL 26/11/08

%Initial organic matter present in the profile expressed in M L-2
%C_manure N_manure C_litter N_litter C_humus N_humus



% soil_om= ones(ncs,1)*[200/100000/100,0.1,200/100000/100,0.500,1000/100000/100,0.001];               
soil_om=ones(ncs,1)*[0.001/10^5,0.1,0.001/10^5,0.500,1000/100000/100,0.001];               
soil_om(:,2) =  soil_om(:,1)/8;
soil_om(:,4) =  soil_om(:,3)/8;
soil_om(:,6) =  soil_om(:,5)/8;
soil_om(:,:) =  soil_om(:,:)*10; %Correction Joachim to avoid extreme (used to be 3)
soil_om(40:60,5:6) =soil_om(40:60,5:6)/10;
soil_om(61:85,5:6) =soil_om(61:85,5:6)/100;
% mineralisaton
 %soil_om= ones(ncs,1)*[1,0.1,5,0.500,0.1,0.001];               
 %soil_om(:,:) = 0;
 soil_om(:,7) = 0; %CO2 (Added by Joachim to have a carbon balance)

 

 
 
 
 