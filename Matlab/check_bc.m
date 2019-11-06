function [top_OK, bot_OK,bctop_changed,bcbot_changed,phsurf,phbot,flxsbot,boco_top_type,...
      boco_bot_type,boco_top]=check_bc(wat_flxs,ph,kh_in,ncs,dx_inter,pond,phsa,...
   	flxar,boco_top_type,boco_bot_type,seep,phsurf,phbot,flxsbot,bctop_changed,...
      bcbot_changed,boco_top); 
   
%Check if BC changed during iteration
%
%IN:
%	wat_flxs,ph,kh_in,compartiments_number,
%	dx_inter,pond,phsa,boco_top_type,boco_top,seep
%   flxar:prescribed upper flux
%OUT:
%	top_OK, bot_OK,bctop_changed
%	phsurf phbot flxsbot
%CALL:
%	none
%CALLED BY:
%   solveflow.m
%--------------------
%M. Javaux 15/05/00

top_OK=1;

if boco_top_type==2
   %max. infiltration
   flxmin = -kh_in(1)*((pond-ph(1))/dx_inter(1)+1);
   %max. evaporation
   flxmax = -kh_in(1)*((phsa-ph(1))/dx_inter(1)+1);
end
if boco_top_type==1 & phsurf<0 & (wat_flxs(1)>flxar)
   boco_top_type=2;
   boco_top = flxar;
   %disp('changed UBC (1)');
   bctop_changed=1;
   top_OK=0;
elseif boco_top_type==1 & phsurf>0 & (wat_flxs(1)<flxar)
   boco_top_type=2;		%infiltration greater than available water => ponding disappears
   top_OK=0;
   bctop_changed=1;
   %disp('changed UBC (2)');
elseif boco_top_type==2 & (flxar>flxmax)	
   boco_top_type=1; % real evaporation greater than max evaporation
   phsurf=phsa;
   top_OK=0;
   bctop_changed=1;
   boco_top = phsa;
   %disp('changed UBC (3)');
elseif boco_top_type==2 & (flxar<flxmin)
   boco_top_type=1;
   phsurf=pond; 			% ponding appears
   %disp('changed UBC (4)')
   bctop_changed=1;
   top_OK=0;
end

% Check bottom boundary condition
bot_OK=1;
if seep==1				%lysimeter bottom boundary
   if boco_bot_type==1 & wat_flxs(ncs+1)>0
      boco_bot_type=2;
      flxsbot=0;
      bot_OK=0;
      bcbot_changed=1;
   elseif boco_bot_type==2 & ph(ncs)>5
      boco_bot_type=1;
      phbot=0;
      bot_OK=0;
      bcbot_changed=1;
   end
end