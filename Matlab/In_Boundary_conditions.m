function [soilConsParams,soilCilinderStateParams]=In_Boundary_conditions(BC_input_option,climateState,soilConsParams,soilCilinderStateParams)
%% DOCUMENTATION
% Inputfile for water boundary conditions
% Enter your top and bottom boundary conditions in a matrix called:
%   BOUNDARY_CONDITIONS_MATRIX
% This matrix is composed by:
% Col 1: time (T)
% Col 2: type of top boundary condition
%		1 = pressure head condition (L)
%		2 = flow condition (L/T): + for evaporation, - for drainage and/or irrigation
% Col 3: top boundary condition
% Col 4: type of bottom boundary condition
%		1 = pressure head condition (L)
%		2 = flow condition (L/T): + for capillarity rise, - for drainage 
%		3 = seepage face (lysimeter)
%		4 = free drainage
%%%AJOUT MAMADOU
%		5 = a groundwater is present and the groundwater level (gwl) is given at each time
%		6 = a groundwater is present and the fonction between gwl and flux is available
%		7 = a groundwater is present and the flux is given at each time
%
% Col 5: bottom boundary condition
% 		if Col 4 = 3, then Col 5 must be equal to 0
%       if Col 4 = 5, then Col 5 = gwl (L) in absolute value
%       if Col 4 = 6, then Col 5 = initial gwl (L)
%       if Col 4 = 7, then Col 5 = flux (L T-1)
%%%%%%%%%%%%%%%%%%%%%%
% IN :
%	none
% OUT :
% 	BOUNDARY_CONDITIONS_MATRIX
%	hs:maximum ponding depth (L)
%	phsa:Minimum allowed pressure head at the surface (L)
%CALLS:
%	none
%
%%%%%%%%%%%%%%%%%%%%%%%
%BC_input_option
% 0 : Raw input of the matrix
% 1 : Use a climatic file and choose appropiate bottom BC
% 2 : Use climatic file for rain and irrigation. ET comes from model 

%% FUNCTION INPUT
% soil common parameters
wat_bottom_BC_type = soilConsParams.wat_bottom_BC_type; % type of bottom boundary condition for water flow
%climate data
climate = climateState.climateDaily;

%% FUNCTION MAIN BODY

if BC_input_option == 0
    %Creata a matrix containing the boundary conditions
    BOUNDARY_CONDITIONS_MATRIX = repmat([0 2 -5 3 0],40,1);
    BOUNDARY_CONDITIONS_MATRIX(:,1) = 0:1:39;
elseif BC_input_option ==1
    %Then fill in the climatic data in in_climatic data and choose the
    %bottom boundary condition, enter tabulated values to split E and T in
    %the files in_crop_number
    BOUNDARY_CONDITIONS_MATRIX(:,1) = climate(:,1);
    BOUNDARY_CONDITIONS_MATRIX(:,2) = 2;
    BOUNDARY_CONDITIONS_MATRIX(:,3) =(climate(:,2)-climate(:,3)-climate(:,4));
    BOUNDARY_CONDITIONS_MATRIX(:,4) = wat_bottom_BC_type(1);
    BOUNDARY_CONDITIONS_MATRIX(:,5) = wat_bottom_BC_type(2);
    
    % FUNCTION OUTPUT
    soilCilinderStateParams.IPBCM = BOUNDARY_CONDITIONS_MATRIX;
    
else
    BOUNDARY_CONDITIONS_MATRIX(:,1) = climate(:,1);
    BOUNDARY_CONDITIONS_MATRIX(:,2) = 2;
    BOUNDARY_CONDITIONS_MATRIX(:,3) = -climate(:,3)-climate(:,4); %Only rain and irrigation, ET comes from crop model
    BOUNDARY_CONDITIONS_MATRIX(:,4) = wat_bottom_BC_type(1);
    BOUNDARY_CONDITIONS_MATRIX(:,5) = wat_bottom_BC_type(2);
      
    % FUNCTION OUTPUT
    soilCilinderStateParams.BOUNDARY_CONDITIONS_MATRIX = BOUNDARY_CONDITIONS_MATRIX;
end

hs=0;								% Maximum ponding depth (L), above this pressure head runoff will occur
phsa=-10000; %origine				% Minimum allowed pressure head at the surface (L)

%Specify the initial gwl
initial_gwl=95;
% for a groundwater flux relationship as bottom condition give the
% coefficient used in the relation Q = arel*exp(brel*gwl)
arel= -0.01;
brel= -0.005;

%% FUNCTION OUTPUT
%soil constant parameters
soilConsParams.arel = arel; %coefficient for groundwater flux
soilConsParams.brel = brel; %coefficient for groundwater flux 
soilConsParams.hs = hs; %maximum pressure head before runoff
soilConsParams.initial_gwl = initial_gwl; % initial groundwater level
soilConsParams.phsa = phsa; % minimum allowed pressure head at the surface


