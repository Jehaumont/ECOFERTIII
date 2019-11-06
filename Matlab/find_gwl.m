function [ncbot,boco_bot,gwl,BOUNDARY_CONDITIONS_MATRIX]=...
    find_gwl(A,B,ph,dx_inter,ncs,BOUNDARY_CONDITIONS_MATRIX,phsurf,arel,brel)
%% DOCUMENTATION
% find the groundwater level (gwl), the node above gwl and the bottom boundary value
%
% ncbot: first unsatuturated node above gwl 
% gwl: groundwater level
%
%CALLED BY Soil_boundary_condition_new
%CALL None
% -------------------
%M. Sall 25/11/09



if A==5  
    gwl= B;
    dxn=cumsum(dx_inter);
    i=max(find(dxn <gwl));
    ncomp=i;
    if ncomp == ncs
        error('ground water level reached bottom of profile')
    end
    ph(ncomp+1)=-abs(gwl)+abs(dxn(ncomp+1)); %%determine ph(ncomp+1) assuming unit gradient
    ncbot=ncomp;
    boco_bot=ph(ncomp+1);               %% Dirichlet condition
    %bocobot=0;ncbot=94;
end
  
%%%%%%Calculation of gwl in the case of boco_bot_type6 or boco_bot_type7%%%%%(calgwl dans WAVE.FOR)
if A ==6 || A ==7
% locate first unsaturated compartment
ncomp=ncs;
    while ph(ncomp)>=0
        if ncomp >=1
            ncomp = ncomp-1;
        end
    end
    if ncomp == ncs%
		error ('ground water level reached bottom of profile')
    end
%%%%%% if first node is saturated  
   if ncomp == 0
       if phsurf >=0
           gwl=phsurf;
       else
           gwl=dxinter(1)/(ph(1)-phsurf)*(phsurf);
       end
%%%%% else calculate depth of first unsaturated node and distance from that node to the groundwater level   
   else
       dunsat=0;
       for i = 1:ncomp
           dunsat=dunsat+dx_inter(i);
       end
       if (ph(ncomp+1)~=ph(ncomp))
           dist=dx_inter(ncomp+1)/(ph(ncomp+1)-ph(ncomp))*(-ph(ncomp));
       else
           dist = 0;
       end
       gwl=(dunsat+dist); 
   end
end

%qdeep is calculated from the relation with gwl (initial gwl is given the others are calculated)
if A==6              
   %arel= -0.01; brel = -0.005;
   qdeep=arel*exp(brel*abs(gwl));
   ncbot = ncomp;
   boco_bot=qdeep;        %% Neuman condition   
   %BOUNDARY_CONDITIONS_MATRIX(i,5)=gwl
end

%%flxsbot=boco_bot is given at each day (in the scrip %%IN_BOUNDARY_CONDITIONS)
if A==7       
    ncbot = ncomp;
    boco_bot=B;         %% Neuman condition 
end