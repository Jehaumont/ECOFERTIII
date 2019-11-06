function [phsurf,boco_top_type]=fix_uboco(phsa,pond,flxar);

% in case of two unsuccesful trials for top boundary condition: 
% switch to dirichlet conditions and repeat loop without checking
%
% IN:
%   phsa,pond,flxar
% OUT:
%   phsurf,boco_top_type
% CALLS:
%   none
% CALLED BY:
%   solve_flow.m
%-------------
% Javaux 10/00

if flxar>0
   boco_top_type=1;
   phsurf=phsa;
else
   boco_top_type=1;
   phsurf=pond;
	disp('flaxar<0')
end
