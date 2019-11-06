function [wat_flxs,iter] = calc_fluxes (ph,kh_in,dt,phsurf,pond,flxsbot,phbot,flxar,...
   boco_top_type,boco_bot_type,maxiter,dx_inter,compartiments_number,iter)

% Calculates the soil water fluxes across the soil nodes
%
%
%IN: 
%	dt,dx_inter, maxiter, compartiments_number 
%	ph , kh_in,phbot,flxsbot,pond,phsurf,flxar
%	boco_top_type,boco_top,boco_bot_type,boco_bot
%OUT:	
%	wat_flxs: flux <0 downward >0 upward
%   iter
%CALL:
%	none
%CALLED BY:
%   solve_flow.m
%
%----------------------------------
% M. Vanclooster 2/2/2000
% modified by M. Javaux 14/05/00
% modified by M. Sall 25/11/09


ncbot = compartiments_number; 

% Surface flux
if boco_top_type == 1												%  Imposed pressure head condition	
   wat_flxs(1)= -kh_in(1)*((phsurf-ph(1))/dx_inter(1)+1); 			% Darcy equation
else 																% boco_top_type == 2												
   wat_flxs(1)=flxar;
   if flxar<0 & pond<ph(1)
   %infiltration and pond <ph(1): ANORMAL so skip to next step diminution
   %iter=maxiter;
   end
end

% Intermediate flux
 wat_flxs(2:ncbot)=-kh_in(2:ncbot).*(-diff(ph)./dx_inter(2:ncbot)+1);


% Bottom flux
if boco_bot_type  == 1 | boco_bot_type  == 5															% Pressure head condition 
   wat_flxs(ncbot+1)=-kh_in(ncbot+1)*((ph(ncbot)-phbot)/dx_inter(ncbot+1)+1);
end    

if boco_bot_type  == 2	| boco_bot_type == 4 | boco_bot_type  == 6 | boco_bot_type  == 7						% Flux condition 
   wat_flxs(ncbot+1)=flxsbot;
end