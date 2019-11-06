function [conc,solute_applic,solboco_top_type,solboco_input_type]=In_solute_boundary(t,nsol)

%Inputfile for solute boundary conditions

%find the actual flux concentration
%IN:
% t: time information (T) 
% nsol: number of solutes
%OUT
% conc: solute concentration (ML-3)
% solboco_top_type: type of top BC for solute

% CALL:none

%CALLED BY:solve_solute
%--------------------------------
% M. Sall 14/04/09  

% Define the type of top boundary condition
% 1 = First type BC = Concentration BC (possible errors on mass
% balance)
% 2 = Third type BC = Concentration Flux BC (no mass balance error)
solboco_top_type = 2; %JV was 2

%Define the type of input : M L-3 or M L-2
% 1 = M L-3
% 2 = M L-2 ==> define details of application in in_solute_upper.m
solboco_input_type = 2; %JV :was 2

%Solute application = [t solute1 solute2 solute3] expressed in M L^-3 and 
%thus the concentration of the water entering the profile  
%if nitrogen is modelled:
% solute 1 = ureum
% solute 2 = ammonium
% solute 3 = nitrate
 solute_applic=...
     [0 000 0 0  
     1 0 0 0
     2  0 0 0]; 
 

   
i=max(find(solute_applic(:,1)<=t));
if numel(i)~=0
    conc=(solute_applic(i,2:end));
else
    conc=zeros(1,nsol);
end
    





