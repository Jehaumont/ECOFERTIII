function [seep,boco_top_type,boco_bot_type,dt,bctop_changed,...
      bcbot_changed,flxar,phbot,phsurf,flxsbot]=calc_boco(dx,kh_in,ph,dt,rtex,t,...
   	compartiments_number,pond,boco_bot,boco_top_type,boco_bot_type,...
      flxa1,stock_max,bctop_changed,bcbot_changed,soil_parameters,flxar,...
      phsa,phsurf,phbot,flxsbot,ponded,pond_from);

%Calculate upper and bottom boundary conditions
%           at every dt step, check if UBC & BCC are right.
%
% IN: 
% 	kh_in:,ph: succion (cm),dt: time step (min)
% 	rtex: root extraction rate in each comp.(1/min), t:time (min)
%	pond: height of ponding (cm), ponded: test if there is ponding
%	flxar: water available at the surface of the soil profile (cm)
%	phsurf: pressure head at the surface of the soil profile (cm)
%	boco_top_type, boco_top, boco_bot_type, boco_bot: type and BC.
% 	bctop_changed, bcbot_changed: test needed by the program "boundary_condition"
% 	dx, stockmax, flxa1
% OUT:
%	seep: test if seepage, phbot: ph of the bottom (cm), phsurf: ph of the surface (cm)
%	flxsbot: flux through the bottom (cm/min),flxar= available flux through the surface (cm/min)
%	boco_top_type, boco_top, boco_bot_type, boco_bot: type and BC.
% 	bctop_changed, bcbot_changed: test needed by the program "boundary_condition"
% CALL:
% 	calc_stock.m
% CALLED BY:
%       solveflow.m
%
%---------------------------
% M. Javaux 15/05/00
%modified by M. Sall 25/11/09
% Bottom boundary condition
%%%%%%%%%%%%%%%%%
ncs=compartiments_number;


if boco_bot_type==1 || boco_bot_type==5						            %ph
    phbot=boco_bot;
    seep=0;
elseif boco_bot_type==2	|| boco_bot_type==6 || boco_bot_type==7				%flux
    flxsbot=boco_bot;
    seep=0;
elseif boco_bot_type==3					%lysimeter
    seep=1;
    if ph(compartiments_number)>=5	%50 comes from Wave.f
        boco_bot_type=1;
        phbot=0;
    else 										%(ph(compartiments_number)<50))
        boco_bot_type=2;
        flxsbot=0;
    end
elseif boco_bot_type==4					% free drainage
    flxsbot=-kh_in(compartiments_number+1);
    seep=0;
end

% Upper boundary condition
% Ponding
if ponded==0 
   if  pond>0
      ponded=1;
      pond_from=t;
   end
end

flxar = flxa1-pond/dt;

% Upper boundary condition
if boco_top_type==1 && phsurf<0						
   if flxar>0
         % flxar >0 => evaporation
      phsurf=phsa;											
   else
      % flxar <0 => infiltration
      boco_top_type=2;
      bctop_changed=1;
      disp('change BC--->flux');
   end
elseif boco_top_type==1 &&  phsurf>=0 
   if flxar>0
      % ponding but evaporation flux
      boco_top_type=2;%define CL evaporation
      bctop_changed=1;
      disp('change BC--->flux');
   else
      %ponding +infiltration
      phsurf=pond;
   end
end

