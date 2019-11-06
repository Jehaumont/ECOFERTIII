function[om_applic,om_appl]=om_app(t,dt,om_appl)
%Dates of application of C and N 
%IN:
%t=time, dt=time increment
%OUT:
%om_appl: amount of organic carbon (C) and organic nitrogen (N) application at time t

%CALL:none
%CALLED BY: miner_immob
%-------------------------------------------------------------------------
% Mamadou SALL 17/03/09
%Correction by Joachim Vansteenkiste 16/08/2010, set application to zero
%when it is applied, avoid double addition in this way.

% time(column1) vs. application of carbon (column2(M L^-2)) and nitrogen(column3(M L-2))
%In miner_immob it is assumed that 30% of the carbon goes to the manure pool, 
%while the rest goes to the litter pool. Manure N addition is derived from
%Manure C addition and ro, Litter N addition is derived from the input
%minus what has already gone in the manure pool.

% if t>163
%     disp('edf')
% end
 JJ=find(floor(t-dt)<=om_appl(:,1) & om_appl(:,1)<=ceil(t));
 if isempty(JJ);
     N_application=0;C_application=0;
 elseif sum(om_appl(JJ,2:3)) == 0;
     N_application=0;C_application=0;     
 else
     C_application=om_appl(JJ,2);
     N_application=om_appl(JJ,3);
     om_appl(JJ,2:3)=0;
 end
 om_applic=[C_application, N_application];
 
 %% Remarks: In Wave.FOR, C_application = carborg and N_application = rnitorg
