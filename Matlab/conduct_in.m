function kh_in=conduct_in(ph,boco_top_type,phsurf,phsa,flxar,soil_parameters)

%Calculation of the hydraulic conductivity in between the nodes
%
%IN:
% 	ph = soil water pressure head (cm)
%   pond,phsa,soil_parameters
%   boco_top_type boco_top boco_bot_type boco_bot
%OUT:
%   kh_in = hydraulic conductivity in between the nodes (m/s)
%CALL:
%	conduct.m
%CALLED BY:
%   solveflow.m
%--------------------------
% M. Vanclooster 2/2/2000

n = length(ph);
kh = conduct(ph,soil_parameters);

if boco_top_type == 1		
   % pressure head boundary condition
   phsurf_est = phsurf ;
end
if boco_top_type == 2			
   % flux type boundary condition 
   if flxar  < 0 		
      % infiltration
      phsurf_est = phsurf;
   else 							
      %evaporation 					
      phsurf_est = phsa;
   end
end

KTOP=conduct(phsurf_est,soil_parameters);				
%=Ksat (infilt) or =ph(1) if ph B.C.

% hydraulic conductivity at the surface
%kh_in(1)=sqrt(KTOP(1)*kh(1));
%aritmecthic mean
kh_in(1)=(KTOP(1)+kh(1))/2;

% hydraulic conductivity inside the soil
 %kh_in(2:n)=sqrt(kh(1:n-1).*kh(2:n));
 %aritmecthic mean
 kh_in(2:n)=(kh(1:n-1)+kh(2:n))/2;
%hydraulic conductivity at the bottom
kh_in(n+1)=kh(n);
